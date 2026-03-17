#########################################
# Title: Label generic Fusions
# Author: Ryley Dorney
# Date: Mar 2026
#########################################
library(AnnotationFilter)
library(ensembldb)
library(GenomicRanges)
library(dplyr)
library(BiocParallel)
library(memoise)

ah <- AnnotationHub()
ensembldbv110 <- ah[["AH113665"]]
#Make read-through/cis-splice check function
check_readthrough <- function(g_left, g_right, edb, id_type = "symbol", chr_left = NULL, chr_right = NULL) {
  # 1. check if input is gene symbol or id, and set-up standard chromosomes
  id_filter <- if(id_type == "symbol") {
    SymbolFilter(c(g_left, g_right)) }
  else {
    GeneIdFilter(c(g_left, g_right))}
  
  std_chrs <- c(1:22, "X", "Y", "MT", "M")
  
  #get the ranges of both genes
  g <- genes(edb, filter = list(id_filter, SeqNameFilter(std_chrs)))
  #ensure both genes are found
  if (length(g) < 2) return(FALSE)
  
  # 2. Get coordinates for both genes separately to preserve order
  
  if(!is.null(chr_left) && !is.null(chr_right)){
    if(id_type == "symbol") {
      left_gene_coord <- genes(edb, filter = list(
        SymbolFilter(g_left),
        SeqNameFilter(chr_left)
      ))
      right_gene_coord <- genes(edb, filter = list(
        SymbolFilter(g_right),
        SeqNameFilter(chr_right)
      ))
    }
    else{ #ID used
      left_gene_coord <- genes(edb, filter = list(
        GeneIdFilter(g_left),
        SeqNameFilter(chr_left)
      ))
      right_gene_coord <- genes(edb, filter = list(
        GeneIdFilter(g_right),
        SeqNameFilter(chr_right)
      ))
    }
  }
  else {#Alternatively, if chromosome is not provided
    if(id_type == "symbol") {
      left_gene_coord  <- g[g$symbol == g_left]
      right_gene_coord <- g[g$symbol == g_right]
    } else {
      left_gene_coord  <- g[g$gene_id == g_left]
      right_gene_coord <- g[g$gene_id == g_right]
    }
  }
  
  # Validation: Must find exactly 2 unique genomic locations for both genes (start and end)
  if (length(left_gene_coord) == 0 || length(right_gene_coord) == 0) return(FALSE)
  if (length(left_gene_coord) > 1 || length(right_gene_coord) > 1) {
    if(!is.null(chr_left) && !is.null(chr_right)){
      # We have multiple hits. Let's force the one that matches our input 'chr_left' and 'chr_right'
      left_gene_coord <- left_gene_coord[seqnames(left_gene_coord) == chr_left]
      right_gene_coord <- right_gene_coord[seqnames(right_gene_coord) == chr_right]
    }
    # Final safety: If there are STILL multiple hits (e.g. same gene, same chr, different IDs)
    # then we take the overall range to avoid the 'length > 1' error.
    if (length(left_gene_coord) > 1 || length(right_gene_coord) > 1) {
      left_gene_coord  <- range(left_gene_coord)
      right_gene_coord <- range(right_gene_coord)
    }
  }
  # 3. Strand and Chromosome Check
  if (as.character(strand(left_gene_coord)[1]) != as.character(strand(right_gene_coord)[1])) return(FALSE)
  if (as.character(seqnames(left_gene_coord)[1]) != as.character(seqnames(right_gene_coord)[1])) return(FALSE)
  
  # 4. Transcriptional Order Check
  # For '+' strand: Left end < Right start
  # For '-' strand: Left start > Right end (coords are higher for 5' gene)
  if (as.character(strand(left_gene_coord)[1]) == "+") {
    # 5' gene end must be less than 3' gene start
    if (end(left_gene_coord) > start(right_gene_coord)) return(FALSE) 
  } else { # check '-' strand
    # 5' gene start (lower num on - strand) must be greater than 3' gene end
    if (start(left_gene_coord) < end(right_gene_coord)) return(FALSE) 
  }
  #5. Intervening Gene Check (NAMED Protein-Coding Only)
  
  combined_range <- range(c(left_gene_coord, right_gene_coord))
  
  if (as.character(strand(left_gene_coord)[1]) == "+") {
    gap_start <- end(left_gene_coord) 
    gap_end <- start(right_gene_coord)
  } else { #'-' strand
    gap_start<- end(right_gene_coord)
    gap_end <- start(left_gene_coord) 
    
  }
  
  gap_range <- GRanges(
    seqnames = seqnames(left_gene_coord),
    ranges   = IRanges(start = gap_start, end = gap_end),
    strand   = strand(left_gene_coord)
  )
  
  intervening <- genes(edb, filter = list(
    GRangesFilter(gap_range),
    SeqStrandFilter(as.character(strand(left_gene_coord)[1])) #,GeneBiotypeFilter("protein_coding")
  ))
  # Remove the query genes themselves 
  intervening <- intervening[!(intervening$symbol %in% c(g_left, g_right))]
  
  #remove genes that overlap the query genes
  query_gr <- c(left_gene_coord, right_gene_coord)
  hits <- findOverlaps(intervening, query_gr, ignore.strand = FALSE)
  
  intervening_clean <- intervening[-queryHits(hits)]
  
  #return(intervening)
  #return(intervening_clean)
  return(length(intervening_clean) == 0)
}
check_SAGe <- function(g_left, g_right, edb, id_type = "symbol", chr_left = NULL, chr_right = NULL) {
  std_chrs <- c(1:22, "X", "Y", "MT", "M")
  
  # 1. check if input is gene symbol or id, and set-up standard chromosomes
  id_filter <- if(id_type == "symbol") {
    SymbolFilter(c(g_left, g_right)) }
  else {
    GeneIdFilter(c(g_left, g_right))}
  
  #get the ranges of both genes
  g <- genes(edb, filter = list(id_filter, SeqNameFilter(std_chrs)))
  #ensure both genes are found
  if (length(g) < 2) return(FALSE)
  
  # 2. Get coordinates for both genes separately to preserve order
  if(!is.null(chr_left) && !is.null(chr_right)){
    if(id_type == "symbol") {
      left_gene_coord <- genes(edb, filter = list(
        SymbolFilter(g_left),
        SeqNameFilter(chr_left)
      ))
      right_gene_coord <- genes(edb, filter = list(
        SymbolFilter(g_right),
        SeqNameFilter(chr_right)
      ))
    }
    else{ #ID used
      left_gene_coord <- genes(edb, filter = list(
        GeneIdFilter(g_left),
        SeqNameFilter(chr_left)
      ))
      right_gene_coord <- genes(edb, filter = list(
        GeneIdFilter(g_right),
        SeqNameFilter(chr_right)
      ))
    }
  }
  else {#Alternatively, if chromosome is not provided
    if(id_type == "symbol") {
      left_gene_coord  <- g[g$symbol == g_left]
      right_gene_coord <- g[g$symbol == g_right]
    } else {
      left_gene_coord  <- g[g$gene_id == g_left]
      right_gene_coord <- g[g$gene_id == g_right]
    }
  }
  # Validation: Must find exactly 2 unique genomic locations for both genes (start and end)
  if (length(left_gene_coord) == 0 || length(right_gene_coord) == 0) return(FALSE)
  if (length(left_gene_coord) > 1 || length(right_gene_coord) > 1) {
    if(!is.null(chr_left) && !is.null(chr_right)){
      # We have multiple hits. Let's force the one that matches our input 'chr_left' and 'chr_right'
      left_gene_coord <- left_gene_coord[seqnames(left_gene_coord) == chr_left]
      right_gene_coord <- right_gene_coord[seqnames(right_gene_coord) == chr_right]
    }
    # Final safety: If there are STILL multiple hits (e.g. same gene, same chr, different IDs)
    # then we take the overall range to avoid the 'length > 1' error.
    if (length(left_gene_coord) > 1 || length(right_gene_coord) > 1) {
      left_gene_coord  <- range(left_gene_coord)
      right_gene_coord <- range(right_gene_coord)
    }
  }
  # 3. Strand and Chromosome Check
  if (as.character(strand(left_gene_coord)[1]) == as.character(strand(right_gene_coord)[1])) return(FALSE)
  if (as.character(seqnames(left_gene_coord)[1]) != as.character(seqnames(right_gene_coord)[1])) return(FALSE)
  
  # 4. Transcriptional Order Check
  # ensure they do not overlap
  if (end(left_gene_coord) >= start(right_gene_coord) && start(left_gene_coord) <= end(right_gene_coord)) {
    return(FALSE)
  }
  
  # define genomic gap (strand-agnostic)
  gap_start <- min(end(left_gene_coord), end(right_gene_coord))
  gap_end   <- max(start(left_gene_coord), start(right_gene_coord))
  
  # if genes touch, gap width = 0 → still adjacent
  if (gap_end > gap_start) {
    
    gap_range <- GRanges(
      seqnames = seqnames(left_gene_coord),
      ranges   = IRanges(start = gap_start, end = gap_end),
      strand   = "*"
    )
    
    # get ANY genes in the gap (both strands)
    intervening <- genes(edb, filter = list(
      GRangesFilter(gap_range)
    ))
    # Remove the query genes themselves 
    intervening <- intervening[!(intervening$symbol %in% c(g_left, g_right))]
    
    #remove genes that overlap the query genes
    query_gr <- c(left_gene_coord, right_gene_coord)
    intervening_clean <- intervening[!intervening %over% query_gr]
    
    #return(intervening)
    #return(intervening_clean)
    return(length(intervening_clean) == 0)
  }
}
m_check_readthrough <- memoise(check_readthrough)
m_check_SAGe        <- memoise(check_SAGe)

#######################
# Label CTAT-LR-Fusion 
#######################
Annot_CTATLR_Huh7 <- CTATLR_sensecheck_annotated

idx <- which(is.na(Annot_CTATLR_Huh7$fusionType) | Annot_CTATLR_Huh7$fusionType == "")# Check if the current fusionType is empty

inter_idx <- idx[Annot_CTATLR_Huh7$chrom1[idx] != Annot_CTATLR_Huh7$chrom2[idx]] #label inter-chromosomal fusions
Annot_CTATLR_Huh7$fusionType[inter_idx] <- "inter-chromosomal"
intra_idx <- idx[Annot_CTATLR_Huh7$chrom1[idx] == Annot_CTATLR_Huh7$chrom2[idx]]

Annot_CTATLR_Huh7$fusionType[intra_idx] <- mapply(function(chr1, chr2, strand1, strand2, gene_name1, gene_name2) {
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    
    if (strand1 == strand2) {
      is_neighbour <- m_check_readthrough(gene_name1, gene_name2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
    } else if (strand1 != strand2){
      is_SAGe <- m_check_SAGe(gene_name1, gene_name2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
    } 
    return("intra-chromosomal")
  
},  
Annot_CTATLR_Huh7$chrom1[intra_idx] , Annot_CTATLR_Huh7$chrom2[intra_idx] ,
Annot_CTATLR_Huh7$strand1[intra_idx] , Annot_CTATLR_Huh7$strand2[intra_idx] , 
Annot_CTATLR_Huh7$LeftGene[intra_idx] , Annot_CTATLR_Huh7$RightGene[intra_idx])

write_tsv(Annot_CTATLR_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/CTATLR_Huh7.tsv")
#######################
# Label Genion 
#######################
# Setup the parallel parameters (5 workers as requested)
param <- SnowParam(workers = 5, progressbar = TRUE)
db_path <- dbfile(dbconn(ensembldbv110))

Annot_Genion_Huh7 <- Genion_sensecheck_annotated
Annot_Genion_Huh7$fusionType <- bpmapply(function(g1, g2, g3, 
                                                 current_type, 
                                                 chr1, chr2, chr3, chr4,
                                                 gene_name1, gene_name2,
                                                 path_to_db, func1, func2) {

  # Check if worker_db already exists in this worker's environment
  if (!exists("worker_db", envir = .GlobalEnv)) {
    libs <- c("stats","ensembldb", "GenomicRanges", "dplyr", "AnnotationFilter")
    lapply(libs, library, character.only = TRUE) 
    
    assign("worker_db", EnsDb(path_to_db), envir = .GlobalEnv)
  }
  # Now use get("worker_db", envir = .GlobalEnv) 
  edb <- get("worker_db", envir = .GlobalEnv)
  
  
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for multi-fusion
    if (!is.na(chr4)) return("tetra-fusion")
    #check for tri-fusion
    if (!is.na(g3) && is.na(chr4)) return("tri-fusion")
    
    # Check for inter-chromosomal first
    if (chr1 != chr2) return("inter-chromosomal") 
    
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    
    if (chr1 == chr2) {
      chromosome <-  chr1   
      is_neighbour <- func1(g1, g2, edb, id_type = "geneid", chr_left = chromosome, chr_right = chromosome)
      is_SAGe <- func2(g1, g2, edb, id_type = "geneid", chr_left = chromosome, chr_right = chromosome) 
      if (is_neighbour) return("read-through")
      if (is_SAGe) return("SAGe") 
      
      return("intra-chromosomal")
    } 
  }
  return(current_type)
} , 
Annot_Genion_Huh7$V1.1, Annot_Genion_Huh7$V1.2,  Annot_Genion_Huh7$V1.3, 
Annot_Genion_Huh7$fusionType, 
Annot_Genion_Huh7$chr1, Annot_Genion_Huh7$chr2 , Annot_Genion_Huh7$chr3, Annot_Genion_Huh7$chr4, 
Annot_Genion_Huh7$V2.1, Annot_Genion_Huh7$V2.2,
MoreArgs = list(path_to_db = db_path, 
                func1 = check_readthrough, # Pass actual function objects
                func2 = check_SAGe),
BPPARAM = param)

write_tsv(Annot_Genion_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/Genion_Huh7.tsv")

# 3. Clean up
bpstop(param) 

#######################
# Label LongGF 
#######################
# Setup the parallel parameters (5 workers as requested)
param <- SnowParam(workers = 5, progressbar = TRUE)
db_path <- dbfile(dbconn(ensembldbv110))

Annot_LongGF_Huh7 <- LongGF_sensecheck_annotated

idx <- which(is.na(Annot_LongGF_Huh7$fusionType) | Annot_LongGF_Huh7$fusionType == "")# Check if the current fusionType is empty
inter_idx <- idx[Annot_LongGF_Huh7$chromosome1[idx] != Annot_LongGF_Huh7$chromosome2[idx]] #label inter-chromosomal fusions
Annot_LongGF_Huh7$fusionType[inter_idx] <- "inter-chromosomal"
intra_idx <- idx[Annot_LongGF_Huh7$chromosome1[idx] == Annot_LongGF_Huh7$chromosome2[idx]]
length(Annot_LongGF_Huh7$fusionType[intra_idx])
length(unique(paste(Annot_LongGF_Huh7$Gene1[intra_idx], Annot_LongGF_Huh7$Gene2[intra_idx])))

Annot_LongGF_Huh7$fusionType[intra_idx] <-  bpmapply(function(g1, g2, chr1, chr2,
                                                              path_to_db, func1, func2) {
  # Check if worker_db already exists in this worker's environment
  if (!exists("worker_db", envir = .GlobalEnv)) {
    libs <- c("stats","ensembldb", "GenomicRanges", "dplyr", "AnnotationFilter")
    lapply(libs, library, character.only = TRUE) 
    
    assign("worker_db", EnsDb(path_to_db), envir = .GlobalEnv)
  }
  # Now use get("worker_db", envir = .GlobalEnv) 
  edb <- get("worker_db", envir = .GlobalEnv)
  
  # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
      is_neighbour <- func1(g1, g2, edb, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- func2(g1, g2, edb, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
  }
 , Annot_LongGF_Huh7$Gene1[intra_idx], Annot_LongGF_Huh7$Gene2[intra_idx], 
Annot_LongGF_Huh7$chromosome1[intra_idx], Annot_LongGF_Huh7$chromosome2[intra_idx],
MoreArgs = list(path_to_db = db_path, 
                func1 = check_readthrough, # Pass actual function objects
                func2 = check_SAGe),
BPPARAM = param)

write_tsv(Annot_LongGF_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/LongGF_Huh7.tsv")

bpstop(param)
#######################
# Label FusionSeeker 
#######################
# Setup the parallel parameters (5 workers as requested)
param <- SnowParam(workers = 5, progressbar = TRUE)
db_path <- dbfile(dbconn(ensembldbv110))

Annot_FusionSeeker_Huh7 <- FusionSeeker_sensecheck_annotated
idx <- which(is.na(Annot_FusionSeeker_Huh7$fusionType) | Annot_FusionSeeker_Huh7$fusionType == "")

inter_idx <- idx[Annot_FusionSeeker_Huh7$Chrom1[idx] != Annot_FusionSeeker_Huh7$Chrom2[idx]] #label inter-chromosomal fusions
Annot_FusionSeeker_Huh7$fusionType[inter_idx] <- "inter-chromosomal"
length(Annot_FusionSeeker_Huh7$fusionType[inter_idx])
intra_idx <- idx[Annot_FusionSeeker_Huh7$Chrom1[idx] == Annot_FusionSeeker_Huh7$Chrom2[idx]]
length(Annot_FusionSeeker_Huh7$fusionType[intra_idx])
length(unique(paste(Annot_FusionSeeker_Huh7$Gene1[intra_idx], Annot_FusionSeeker_Huh7$Gene2[intra_idx])))

Annot_FusionSeeker_Huh7$fusionType[intra_idx] <- bpmapply(
  function(g1, g2, chr1, chr2,
           path_to_db, func1, func2) {
    # Check if worker_db already exists in this worker's environment
    if (!exists("worker_db", envir = .GlobalEnv)) {
      libs <- c("stats","ensembldb", "GenomicRanges", "dplyr", "AnnotationFilter")
      lapply(libs, library, character.only = TRUE) 
      
      assign("worker_db", EnsDb(path_to_db), envir = .GlobalEnv)
    }
    # Now use get("worker_db", envir = .GlobalEnv) 
    edb <- get("worker_db", envir = .GlobalEnv)
    
    #check types
    chromosome <- sub("^chr", "", chr1, ignore.case = TRUE)
    
    is_neighbour <- func1(
      g1, g2, edb,
      id_type="geneid",
      chr_left=chromosome,
      chr_right=chromosome
    )
    if (is_neighbour) return("read-through")
    
    is_SAGe <- func2(
      g1, g2, edb,
      id_type="geneid",
      chr_left=chromosome,
      chr_right=chromosome
    )
    if (is_SAGe) return("SAGe")
    
    "intra-chromosomal"
    
  },
  Annot_FusionSeeker_Huh7$Gene1[intra_idx],
  Annot_FusionSeeker_Huh7$Gene2[intra_idx],
  Annot_FusionSeeker_Huh7$Chrom1[intra_idx],
  Annot_FusionSeeker_Huh7$Chrom2[intra_idx],
  MoreArgs = list(path_to_db = db_path, 
                  func1 = m_check_readthrough, # Pass actual function objects
                  func2 = m_check_SAGe),
  BPPARAM = param
)

write_tsv(Annot_FusionSeeker_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/FusionSeeker_Huh7.tsv")
bpstop(param)

#######################
# Label GFSeeker 
#######################
# Setup the parallel parameters (5 workers as requested)
param <- SnowParam(workers = 5, progressbar = TRUE)
db_path <- dbfile(dbconn(ensembldbv110))

Annot_GFSeeker_Huh7 <- GFSeeker_sensecheck_annotated

idx <- which(is.na(Annot_GFSeeker_Huh7$fusionType) | Annot_GFSeeker_Huh7$fusionType == "")

inter_idx <- idx[Annot_GFSeeker_Huh7$chrom1[idx] != Annot_GFSeeker_Huh7$chrom2[idx]] #label inter-chromosomal fusions
Annot_GFSeeker_Huh7$fusionType[inter_idx] <- "inter-chromosomal"
length(Annot_GFSeeker_Huh7$fusionType[inter_idx])
intra_idx <- idx[Annot_GFSeeker_Huh7$chrom1[idx] == Annot_GFSeeker_Huh7$chrom2[idx]]
length(Annot_GFSeeker_Huh7$fusionType[intra_idx])
length(unique(paste(Annot_GFSeeker_Huh7$gene1_name[intra_idx], Annot_GFSeeker_Huh7$gene2_name[intra_idx])))


Annot_GFSeeker_Huh7$fusionType[intra_idx] <- bpmapply(function(g1, g2,chr1, chr2,
                                                               path_to_db, func1, func2) {
  # Check if worker_db already exists in this worker's environment
  if (!exists("worker_db", envir = .GlobalEnv)) {
    libs <- c("stats","ensembldb", "GenomicRanges", "dplyr", "AnnotationFilter")
    lapply(libs, library, character.only = TRUE) 
    
    assign("worker_db", EnsDb(path_to_db), envir = .GlobalEnv)
  }
  # Now use get("worker_db", envir = .GlobalEnv) 
  edb <- get("worker_db", envir = .GlobalEnv)
  
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2){
      is_neighbour <- func1(g1, g2, edb, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- func2(g1, g2, edb, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
    }

} , 
Annot_GFSeeker_Huh7$gene1_name[intra_idx], Annot_GFSeeker_Huh7$gene2_name[intra_idx], 
Annot_GFSeeker_Huh7$chrom1[intra_idx], Annot_GFSeeker_Huh7$chrom2[intra_idx],
MoreArgs = list(path_to_db = db_path, 
                func1 = m_check_readthrough, # Pass actual function objects
                func2 = m_check_SAGe),
BPPARAM = param)

write_tsv(Annot_GFSeeker_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/GFSeeker_Huh7.tsv")

bpstop(param)
#######################
# Label JAFFAL 
#######################
#Repeat for tri-gene fusions
Annot_JAFFAL_Huh7_3Gene <- Huh7_JAFFAL_3Gene
Annot_JAFFAL_Huh7_3Gene$fusionType <- "tri-fusion"
Annot_JAFFAL_Huh7_3Gene$Algorithm <- "JAFFAL"
write_tsv(Annot_JAFFAL_Huh7_3Gene, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/JAFFAL_3Gene_Huh7.tsv")

# Setup the parallel parameters (5 workers as requested)
ensembldbv109 <- ah[["AH109606"]]
param <- SnowParam(workers = 5, progressbar = TRUE)
db_path <- dbfile(dbconn(ensembldbv109))

Annot_JAFFAL_Huh7 <- JAFFAL_sensecheck_annotated

idx <- which(is.na(Annot_JAFFAL_Huh7$fusionType) | Annot_JAFFAL_Huh7$fusionType == "")# Check if the current fusionType is empty
               
inter_idx <- idx[Annot_JAFFAL_Huh7$chrom1[idx] != Annot_JAFFAL_Huh7$chrom2[idx]] #label inter-chromosomal fusions
length(Annot_JAFFAL_Huh7$fusionType[inter_idx])
Annot_JAFFAL_Huh7$fusionType[inter_idx] <- "inter-chromosomal" 

#label intra-chromosomal fusions and their sub
intra_idx <- idx[Annot_JAFFAL_Huh7$chrom1[idx] == Annot_JAFFAL_Huh7$chrom2[idx]]
length(Annot_JAFFAL_Huh7$fusionType[intra_idx])
length(unique(paste(Annot_JAFFAL_Huh7$Gene1[intra_idx], Annot_JAFFAL_Huh7$Gene2[intra_idx])))
length(unique(paste(Annot_JAFFAL_Huh7$ensembl_gene_id.x[intra_idx], Annot_JAFFAL_Huh7$ensembl_gene_id.y[intra_idx])))

Annot_JAFFAL_Huh7$fusionType[intra_idx] <- bpmapply(function(g1, g2, id1, id2, chr1, chr2, str1, str2,
                                                             path_to_db, func1, func2) {
  # Check if worker_db already exists in this worker's environment
  if (!exists("worker_db", envir = .GlobalEnv)) {
    libs <- c("stats","ensembldb", "GenomicRanges", "dplyr", "AnnotationFilter")
    lapply(libs, library, character.only = TRUE) 
    
    assign("worker_db", EnsDb(path_to_db), envir = .GlobalEnv)
  }
  # Now use get("worker_db", envir = .GlobalEnv) 
  edb <- get("worker_db", envir = .GlobalEnv)
  
  # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2) {
      if (str1 == str2) {
      is_neighbour <- func1(id1, id2, edb, id_type = "geneid", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      } else if (str1 != str2) {
      is_SAGe <- func2(id1, id2, edb, id_type = "geneid", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      }
      return("intra-chromosomal")
    } 
} , Annot_JAFFAL_Huh7$Gene1[intra_idx], Annot_JAFFAL_Huh7$Gene2[intra_idx], 
Annot_JAFFAL_Huh7$ensembl_gene_id.x[intra_idx], Annot_JAFFAL_Huh7$ensembl_gene_id.y[intra_idx],
Annot_JAFFAL_Huh7$chrom1[intra_idx], Annot_JAFFAL_Huh7$chrom2[intra_idx],
Annot_JAFFAL_Huh7$strand1[intra_idx], Annot_JAFFAL_Huh7$strand2[intra_idx],
MoreArgs = list(path_to_db = db_path, 
                func1 = m_check_readthrough, # Pass actual function objects
                func2 = m_check_SAGe),
BPPARAM = param) 

write_tsv(Annot_JAFFAL_Huh7, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/JAFFAL_JAFFAdirect_Huh7.tsv")

length(which(Annot_JAFFAL_Huh7$fusionType == "read-through" | Annot_JAFFAL_Huh7$fusionType == "SAGe"))

bpstop(param)
#####################################
# Plot fusion types
#####################################
fusions_Huh7_types <- rbind(dplyr::select(Annot_CTATLR_Huh7, 
                                          c("RNA_sample", "Platform", "Cell_Line", 
                                            "Algorithm", "library_type",    
                                            "fusionType", "full_label")),
                            dplyr::select(Annot_Genion_Huh7, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",   
                                             "fusionType", "full_label")),
                            dplyr::select(Annot_LongGF_Huh7, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",    
                                            "fusionType", "full_label")), 
                            dplyr::select(Annot_FusionSeeker_Huh7, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",   
                                             "fusionType", "full_label")),
                            dplyr::select(Annot_GFSeeker_Huh7, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",    
                                            "fusionType", "full_label")),
                            dplyr::select(Annot_JAFFAL_Huh7, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",    
                                            "fusionType", "full_label")),
                            dplyr::select(Annot_JAFFAL_Huh7_3Gene, 
                                          c("RNA_sample",  "Platform","Cell_Line", 
                                            "Algorithm", "library_type",    
                                            "fusionType", "full_label")))

fusions_Huh7_types$Platform <- factor(fusions_Huh7_types$Platform, levels = c("Illumina", "PacBio", "ONT")) 

fusions_Huh7_types$library_type <- factor(fusions_Huh7_types$library_type, levels = c("PCR_cDNA", "direct_cDNA", "direct_RNA")) 
fusions_Huh7_types$full_label <- factor(fusions_Huh7_types$full_label, levels = c("Illumina.PCR_cDNA.B1", "Illumina.PCR_cDNA.B2", 
                                                                                  "PacBio.PCR_cDNA.B1", "PacBio.PCR_cDNA.B2", 
                                                                                  "ONT.PCR_cDNA.B1","ONT.PCR_cDNA.B2", 
                                                                                  "ONT.direct_cDNA.B1", "ONT.direct_cDNA.B2", 
                                                                                  "ONT.direct_RNA.B1", "ONT.direct_RNA.B2"))

write_tsv(fusions_Huh7_types, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/fusions_Huh7_types.tsv")


ggplot(fusions_Huh7_types, 
       aes(x = RNA_sample, colour = full_label, shape = Algorithm)) +
  geom_point(stat = "count") +
  theme_bw() +
  labs(title = "Fusion Types", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  #labs(fill = "Fusion Type")+
  scale_y_log10()+
  scale_colour_manual(values = platformlibsamp_colourmap)

ggplot(dplyr::filter(fusions_Huh7_types, !fusionType %in% c(
  "Mitochondrial:Genomic",
  "Mitochondrial:Mitochondrial",
  "Self-Misalignment"
)), 
       aes(x = RNA_sample, colour = Algorithm)) +
  geom_point(stat = "count") +
  theme_bw() +
  labs(title = "", subtitle = "minimum read support of 2", x = "", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  labs(fill = "Algorithm")+
  scale_y_log10()+
  scale_colour_manual(values = Alg_colour_map)

library(dplyr)

count_data <- fusions_Huh7_types %>%
  filter(!fusionType %in% c(
    "Mitochondrial:Genomic",
    "Mitochondrial:Mitochondrial",
    "Self-Misalignment"
  )) %>%
  count(RNA_sample, Algorithm, fusionType, library_type, Platform)

ggplot(
  count_data,
  aes(x = RNA_sample, y = n)
) +
  geom_boxplot(aes(group = RNA_sample), outlier.shape = NA) +
  geom_jitter(aes(colour = Algorithm), width = 0.2, size = 2) +
  theme_bw() +
  labs(
    title = "",
    subtitle = "minimum read support of 2",
    x = "",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  facet_grid(fusionType ~ library_type + Platform) +
  scale_y_log10() +
  scale_colour_manual(values = Alg_colour_map)

library(dplyr)
library(tidyr)

paired_data <- count_data %>%
  filter(RNA_sample %in% c("B1", "B2")) %>%
  pivot_wider(
    names_from = RNA_sample,
    values_from = n
  ) %>%
  filter(!is.na(B1) & !is.na(B2))

