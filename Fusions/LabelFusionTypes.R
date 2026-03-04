###############################
# Label other types of fusions
###############################
library(AnnotationFilter)

genes <- gencodev44gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

ah <- AnnotationHub()
ensembldbv110 <- ah[["AH113665"]]
ensembldbv109 <- ah[["AH109606"]]

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
  
  gap_range <- G###############################
# Label other types of fusions
###############################
library(AnnotationFilter)

genes <- gencodev44gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

ah <- AnnotationHub()
ensembldbv110 <- ah[["AH113665"]]
ensembldbv109 <- ah[["AH109606"]]

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
# Example usage:
check_readthrough("CTSD","IFITM10",ensembldbv110, id_type = "symbol", chr_left = "11", chr_right = "11")
check_readthrough("ENSG00000008018","ENSG00000288829",ensembldbv110, id_type = "geneid")
check_readthrough("PTMAP9","MIR1244-4",ensembldbv110, id_type = "symbol")
check_readthrough("ENSG00000198134","ENSG00000283475",ensembldbv110, id_type = "geneid")

#Make read-through/cis-splice check function
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
check_SAGe("CTSD","IFITM10",ensembldbv110, id_type = "symbol")
check_SAGe("CEP70","FAIM",ensembldbv110, id_type = "symbol")

#######################
# Label CTAT-LR-Fusion 
#######################
Annot_CTATLR_Huh7 <- CTATLR_Huh7_sensemito_annotated

Annot_CTATLR_Huh7$fusionType <- mapply(function(current_type, chr1, chr2, strand1, strand2, gene_name1, gene_name2) {
  
  # Process only if current_type is empty or NA
  if (is.na(current_type) || current_type == "") {
    
    # Check for inter-chromosomal first
    if (chr1 != chr2) {
      return("inter-chromosomal") 
    } 
    
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    
    if (strand1 == strand2) {
      is_neighbour <- check_readthrough(gene_name1, gene_name2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
    } else { #(strand1 != strand2)
      is_SAGe <- check_SAGe(gene_name1, gene_name2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
    } 
    return("intra-chromosomal")
  }
  
  # If fusionType was already set (e.g., 'distal-conflicting'), keep it
  return(current_type)
  
},  
Annot_CTATLR_Huh7$fusionType, 
Annot_CTATLR_Huh7$chrom1, Annot_CTATLR_Huh7$chrom2,
Annot_CTATLR_Huh7$strand1, Annot_CTATLR_Huh7$strand2, 
Annot_CTATLR_Huh7$LeftGene, Annot_CTATLR_Huh7$RightGene)


#######################
# Label Genion 
#######################
Annot_Genion_Huh7 <- Genion_Huh7_sensemito_annotated

Annot_Genion_Huh7$fusionType <- mapply(function(g1, g2, g3, 
                                                 current_type, 
                                                 chr1, chr2, chr3, 
                                                 gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for tri-fusion
    if (!is.na(g3)) return("tri-fusion")
    
    # Check for inter-chromosomal first
    if (chr1 != chr2) return("inter-chromosomal") 
    
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    
    if (chr1 == chr2) {
      chromosome <-  chr1   
      is_neighbour <- check_readthrough(g1, g2, ensembldbv110, id_type = "geneid", chr_left = chromosome, chr_right = chromosome)
      is_SAGe <- check_SAGe(g1, g2, ensembldbv110, id_type = "geneid", chr_left = chromosome, chr_right = chromosome) 
      if (is_neighbour) return("read-through")
      if (is_SAGe) return("SAGe") 
      
      return("intra-chromosomal")
    } 
  }
  return(current_type)
} , Annot_Genion_Huh7$V1.1, Annot_Genion_Huh7$V1.2,  Annot_Genion_Huh7$V1.3, 
Annot_Genion_Huh7$fusionType, 
Annot_Genion_Huh7$chr1, Annot_Genion_Huh7$chr2 , Annot_Genion_Huh7$chr3, 
Annot_Genion_Huh7$V2.1, Annot_Genion_Huh7$V2.2)

#######################
# Label LongGF 
#######################
Annot_LongGF_Huh7 <- LongGF_Huh7_sensemito_annotated
Annot_LongGF_Huh7$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2){
      is_neighbour <- check_readthrough(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- check_SAGe(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
    }
  }
  
  # If fusionType was already set (e.g., 'distal-conflicting'), keep it
  return(current_type)
} , Annot_LongGF_Huh7$Gene1, Annot_LongGF_Huh7$Gene2, 
Annot_LongGF_Huh7$fusionType, 
Annot_LongGF_Huh7$chromosome1, Annot_LongGF_Huh7$chromosome2)

#######################
# Label FusionSeeker 
#######################
Annot_FusionSeeker_Huh7 <- FusionSeeker_Huh7_sensemito_annotated
Annot_FusionSeeker_Huh7$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    # Check for inter-chromosomal first
    if (chr1 != chr2) {
      return("inter-chromosomal") 
    } 
    
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2) {
      is_neighbour <- check_readthrough(g1, g2, ensembldbv110, id_type = "geneid", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- check_SAGe(g1, g2, ensembldbv110, id_type = "geneid", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
    } 
  }
  return(current_type)
} , Annot_FusionSeeker_Huh7$Gene1, Annot_FusionSeeker_Huh7$Gene2, 
Annot_FusionSeeker_Huh7$fusionType, 
Annot_FusionSeeker_Huh7$Chrom1, Annot_FusionSeeker_Huh7$Chrom2)

#######################
# Label GFSeeker 
#######################
Annot_GFSeeker_Huh7 <- GFSeeker_Huh7_sensemito_annotated
Annot_GFSeeker_Huh7$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2){
      is_neighbour <- check_readthrough(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- check_SAGe(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
    }
  }
  return(current_type)
} , Annot_GFSeeker_Huh7$gene1_name, Annot_GFSeeker_Huh7$gene2_name, 
Annot_GFSeeker_Huh7$fusionType, 
Annot_GFSeeker_Huh7$chrom1, Annot_GFSeeker_Huh7$chrom2)

#######################
# Label JAFFAL 
#######################
Annot_JAFFAL_Huh7 <- JAFFAL_Huh7_sensemito_annotated
Annot_JAFFAL_Huh7$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    # Check for inter-chromosomal first
    if (chr1 != chr2) {
      return("inter-chromosomal") 
    } 
    
    # Handle intra-chromosomal
    # This calls the function defined above
    # It will return TRUE if they are direct neighbors on the same strand
    chromosome <- sub("^chr", "", chr1,  ignore.case = TRUE)
    if (chr1 == chr2) {
      is_neighbour <- check_readthrough(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome)
      if (is_neighbour) { return("read-through") }
      
      is_SAGe <- check_SAGe(g1, g2, ensembldbv110, id_type = "symbol", chr_left = chromosome, chr_right = chromosome) 
      if (is_SAGe) { return("SAGe") } 
      
      return("intra-chromosomal")
    } 
  }
  return(current_type)
} , Annot_JAFFAL_Huh7$Gene1, Annot_JAFFAL_Huh7$Gene2, 
Annot_JAFFAL_Huh7$fusionType, 
Annot_JAFFAL_Huh7$chrom1, Annot_JAFFAL_Huh7$chrom2)

#Repeat for tri-gene fusions
Annot_JAFFAL_Huh7_3Gene <- JAFFAL_Huh7_3Gene
Annot_JAFFAL_Huh7_3Gene$fusionType <- "tri-fusion"
