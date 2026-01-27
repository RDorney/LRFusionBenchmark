###############################
# Check antisense sense fusions
###############################
#Load necessary library
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
#load relevant gtf file
gencodev44gtf <- import("/bioinformatics/ryley/Gencode44/reference_v44/gencode.v44.annotation.gtf")
############################
# set up dataframes to check
############################
CTATLR_sensecheck <- CTATLR_SGNex
JAFFAL_sensecheck <- JAFFAL_SGNex
Genion_sensecheck <- Genion_SGNex
FusionSeeker_sensecheck <- FusionSeeker_SGNex
LongGF_sensecheck <- LongGF_SGNex

########################################
# Make function to find antisense genes
########################################
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

get_antisense_overlaps <- function(gene, genes_gr, id_type = c("gene_name", "gene_id")) {
  
  id_type <- match.arg(id_type)
  
  if (id_type == "gene_id") {
    # strip version if present (ENSG00000xxxx.y)
    gene <- sub("\\..*$", "", gene)
    gene_ids <- sub("\\..*$", "", genes_gr$gene_id)
    target <- genes_gr[gene_ids == gene]
  } else {
    target <- genes_gr[genes_gr$gene_name == gene]
  }
  
  if (length(target) == 0)
    stop("Gene not found in annotation")
  
  hits <- findOverlaps(target, genes_gr, ignore.strand = TRUE)
  
  query_idx   <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  
  hit_df <- data.frame(
    target_gene   = as.character(target$gene_name[query_idx]),
    target_gene_id = as.character(target$gene_id[query_idx]),
    target_strand = as.character(strand(target)[query_idx]),
    gene_name     = as.character(genes_gr$gene_name[subject_idx]),
    gene_id       = as.character(genes_gr$gene_id[subject_idx]),
    gene_type     = genes_gr$gene_type[subject_idx],
    seqnames      = as.character(seqnames(genes_gr)[subject_idx]),
    start         = start(genes_gr)[subject_idx],
    end           = end(genes_gr)[subject_idx],
    strand        = as.character(strand(genes_gr)[subject_idx]),
    stringsAsFactors = FALSE
  )
  
  hit_df %>%
    dplyr::filter(
      gene_id != target_gene_id,   # remove self
      strand != target_strand      # opposite strand
    ) %>%
    dplyr::distinct(
      gene_name, gene_id, gene_type,
      seqnames, start, end, strand
    ) 
}

#get_antisense_overlaps("MALAT1", genes) #check function works
#get_antisense_overlaps("ENSG00000287557", genes) # check function works

###########################
# Check CTAT-LR
###########################
# get unique genes from both columns
genes_to_check <- unique(c(CTATLR_sensecheck$LeftGene, CTATLR_sensecheck$RightGene))

# run get_antisense_overlaps on each gene
CTATLR_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_name"))

# combine into a single data frame
CTATLR_antisense_df <- bind_rows(CTATLR_antisense_list, .id = "query_gene_index")

# add the original gene names
CTATLR_antisense_df$query_gene <- genes_to_check[as.integer(CTATLR_antisense_df$query_gene_index)]
CTATLR_antisense_df$query_gene_index <- NULL

antisense_pairs <- CTATLR_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

CTATLR_sensecheck_annotated <- CTATLR_sensecheck %>%
  mutate(
    fusionType = if_else(
      # Left = query_gene AND Right = gene_name
      paste(LeftGene, RightGene) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
      paste(RightGene, LeftGene) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial fusions
CTATLR_sensemito_annotated <- CTATLR_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      grepl("chrM:", LeftBreakpoint,  ignore.case = TRUE) &
        grepl("chrM:", RightBreakpoint, ignore.case = TRUE)
      ~ "Mitochondrial:Mitochondrial",
      xor(
        grepl("chrM:", LeftBreakpoint,  ignore.case = TRUE),
        grepl("chrM:", RightBreakpoint, ignore.case = TRUE)
      ) ~ "Mitochondrial:Genomic",
      (LeftGene == RightGene) ~ "Self-Misalignment",
      TRUE ~ fusionType
    )
  )

###########################
# Check Genion
###########################
# get unique genes from both columns
genes_to_check <- unique(c(Genion_sensecheck$V1.1, Genion_sensecheck$V1.2))

# run get_antisense_overlaps on each gene
Genion_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_id"))

# combine into a single data frame
Genion_antisense_df <- bind_rows(Genion_antisense_list, .id = "query_gene_index")

# add the original gene names
Genion_antisense_df$query_gene <- genes_to_check[as.integer(Genion_antisense_df$query_gene_index)]
Genion_antisense_df$query_gene_index <- NULL

Genion_sensecheck_annotated <- Genion_sensecheck %>%
  mutate(
    fusionType = if_else(
      # Left = query_gene AND Right = gene_name
      paste(V1.1, V1.2) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(V1.2, V1.1) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "sense-antisense",
      NA_character_
    )
  )



###########################
# Check FusionSeeker
###########################
# get unique genes from both columns
genes_to_check <- unique(c(FusionSeeker_sensecheck$Gene1, FusionSeeker_sensecheck$Gene2))

# run get_antisense_overlaps on each gene
FusionSeeker_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_id"))

# combine into a single data frame
FusionSeeker_antisense_df <- bind_rows(FusionSeeker_antisense_list, .id = "query_gene_index")

# add the original gene names
FusionSeeker_antisense_df$query_gene <- genes_to_check[as.integer(FusionSeeker_antisense_df$query_gene_index)]
FusionSeeker_antisense_df$query_gene_index <- NULL
###########################
# Check LongGF
###########################
# get unique genes from both columns
genes_to_check <- unique(c(LongGF_sensecheck$Gene1, LongGF_sensecheck$Gene2))

# run get_antisense_overlaps on each gene
LongGF_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_name"))

# combine into a single data frame
LongGF_antisense_df <- bind_rows(LongGF_antisense_list, .id = "query_gene_index")

# add the original gene names
LongGF_antisense_df$query_gene <- genes_to_check[as.integer(LongGF_antisense_df$query_gene_index)]
LongGF_antisense_df$query_gene_index <- NULL

##################################
# Check JAFFAL and update function
##################################
gencodev43gtf <- import("/bioinformatics/ryley/reference_files/gencode.v43.annotation.gtf")

genes <- gencodev43gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

get_antisense_overlaps <- function(gene, genes_gr, id_type = c("gene_name", "gene_id")) {
  
  id_type <- match.arg(id_type)
  
  if (id_type == "gene_id") {
    # strip version if present (ENSG00000xxxx.y)
    gene <- sub("\\..*$", "", gene)
    gene_ids <- sub("\\..*$", "", genes_gr$gene_id)
    target <- genes_gr[gene_ids == gene]
  } else {
    target <- genes_gr[genes_gr$gene_name == gene]
  }
  
  if (length(target) == 0)
    stop("Gene not found in annotation")
  
  hits <- findOverlaps(target, genes_gr, ignore.strand = TRUE)
  
  query_idx   <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  
  hit_df <- data.frame(
    target_gene   = as.character(target$gene_name[query_idx]),
    target_gene_id = as.character(target$gene_id[query_idx]),
    target_strand = as.character(strand(target)[query_idx]),
    gene_name     = as.character(genes_gr$gene_name[subject_idx]),
    gene_id       = as.character(genes_gr$gene_id[subject_idx]),
    gene_type     = genes_gr$gene_type[subject_idx],
    seqnames      = as.character(seqnames(genes_gr)[subject_idx]),
    start         = start(genes_gr)[subject_idx],
    end           = end(genes_gr)[subject_idx],
    strand        = as.character(strand(genes_gr)[subject_idx]),
    stringsAsFactors = FALSE
  )
  
  hit_df %>%
    dplyr::filter(
      gene_id != target_gene_id,   # remove self
      strand != target_strand      # opposite strand
    ) %>%
    dplyr::distinct(
      gene_name, gene_id, gene_type,
      seqnames, start, end, strand
    ) 
}

# get unique genes from both columns
genes_to_check <- unique(c(JAFFAL_sensecheck$Gene1, JAFFAL_sensecheck$Gene2))

# run get_antisense_overlaps on each gene
JAFFAL_antisense_list <- lapply(genes_to_check, function(g) {
  if (grepl("^ENSG", g)) {
    get_antisense_overlaps(g, genes, id_type = "gene_id")
  } else {
    get_antisense_overlaps(g, genes, id_type = "gene_name")
  }
})

# combine into a single data frame
JAFFAL_antisense_df <- bind_rows(JAFFAL_antisense_list, .id = "query_gene_index")
# add the original gene names
JAFFAL_antisense_df$query_gene <- genes_to_check[as.integer(JAFFAL_antisense_df$query_gene_index)]
JAFFAL_antisense_df$query_gene_index <- NULL
