###############################
# Check antisense sense fusions
###############################
#Load necessary library
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
#load relevant gtf file
gencodev44gtf <- import("/bioinformatics/ryley/Gencode44/reference_v44/gencode.v44.annotation.gtf")
# Prep files
CTATLR_sensecheck <- CTATLR_Huh7_ensemblID
JAFFAL_sensecheck <- Huh7_JAFFAL_ensemblID
Genion_sensecheck <- Genion_Huh7
FusionSeeker_sensecheck <- FusionSeeker_Huh7_ensemblID
LongGF_sensecheck <- LongGF_Huh7_ensemblID
GFSeeker_sensecheck <- GFSeeker_Huh7_ensemblID
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
# Annotate with mitochondrial and self fusions
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

antisense_pairs <- Genion_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

Genion_sensecheck_annotated <- Genion_sensecheck %>%
  mutate(
    fusionType = if_else(
      paste(V1.1, V1.2) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(V1.2, V1.1) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial and self fusions
Genion_sensemito_annotated <- Genion_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      grepl("M:", chr1,  ignore.case = TRUE) &
        grepl("M:", chr2, ignore.case = TRUE) 
      ~ "Mitochondrial:Mitochondrial",
      (
        grepl("M:", chr1,  ignore.case = TRUE) +
          grepl("M:", chr2, ignore.case = TRUE) +
          grepl("M:", chr3, ignore.case = TRUE)
      ) == 1 ~ "Mitochondrial:Genomic",
      (V1.1 == V1.2) ~ "Self-Misalignment",
      TRUE ~ fusionType
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

antisense_pairs <- FusionSeeker_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

FusionSeeker_sensecheck_annotated <- FusionSeeker_sensecheck %>%
  mutate(
    fusionType = if_else(
      paste(Gene1, Gene2) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(Gene2, Gene1) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial and self fusions
FusionSeeker_sensemito_annotated <- FusionSeeker_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      Chrom1 == "chrM" & Chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == Chrom1 , "chrM" == Chrom2) ~ "Mitochondrial:Genomic",
      Gene1 == Gene2 ~ "Self-Misalignment",
      TRUE ~ fusionType
    )
  )

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

antisense_pairs <- LongGF_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

LongGF_sensecheck_annotated <- LongGF_sensecheck %>%
  mutate(
    fusionType = if_else(
      paste(Gene1, Gene2) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(Gene2, Gene1) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial and self fusions
LongGF_sensemito_annotated <- LongGF_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      chromosome1 == "chrM" & chromosome2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chromosome1, "chrM" == chromosome2) ~ "Mitochondrial:Genomic",
      Gene1 == Gene2 ~ "Self-Misalignment",
      TRUE ~ fusionType
    )
  )

##################################
# Check GFSeeker
##################################
# get unique genes from both columns
genes_to_check <- unique(c(GFSeeker_sensecheck$gene1_name, GFSeeker_sensecheck$gene2_name))

# run get_antisense_overlaps on each gene
GFSeeker_antisense_list <- lapply(genes_to_check, 
                                  function(g) get_antisense_overlaps(g, genes, id_type = "gene_name"))

# combine into a single data frame
GFSeeker_antisense_df <- bind_rows(GFSeeker_antisense_list, .id = "query_gene_index")

# add the original gene names
GFSeeker_antisense_df$query_gene <- genes_to_check[as.integer(GFSeeker_antisense_df$query_gene_index)]
GFSeeker_antisense_df$query_gene_index <- NULL

antisense_pairs <- GFSeeker_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

GFSeeker_sensecheck_annotated <- GFSeeker_sensecheck %>%
  mutate(
    fusionType = if_else(
      paste(gene1_name, gene2_name) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(gene2_name, gene1_name) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial and self fusions
GFSeeker_sensemito_annotated <- GFSeeker_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      chrom1 == "chrM" & chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chrom1, "chrM" == chrom2) ~ "Mitochondrial:Genomic",
      gene1_name == gene2_name ~ "Self-Misalignment",
      TRUE ~ fusionType
    )
  )

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

antisense_pairs <- JAFFAL_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

JAFFAL_sensecheck_annotated <- JAFFAL_sensecheck %>%
  mutate(
    fusionType = if_else(
      paste(Gene1, Gene2) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
      | 
        paste(Gene2, Gene1) %in% 
        paste(antisense_pairs$query_gene, antisense_pairs$gene_name),
      "Sense-Antisense",
      NA_character_
    )
  )
# Annotate with mitochondrial and self fusions
JAFFAL_sensemito_annotated <- JAFFAL_sensecheck_annotated %>%
  mutate(
    fusionType = case_when(
      chrom1 == "chrM" & chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chrom1, "chrM" == chrom2) ~ "Mitochondrial:Genomic",
      (Gene1 == Gene2) ~ "Self-Misalignment",
      TRUE ~ fusionType
    )
  )

###########################
# combined dataframes
###########################
sensemito_fusions <- rbind(dplyr::select(CTATLR_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")),
                           dplyr::select(Genion_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")),
                           dplyr::select(LongGF_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")), 
                           dplyr::select(FusionSeeker_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")),
                           dplyr::select(GFSeeker_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")),
                           dplyr::select(JAFFAL_sensemito_annotated, 
                                         c("Source", "Cell_Lines", 
                                           "Algorithm", "Sequencing_Depth",    
                                           "Library", "fusionType")))

################
# plots
################
ggplot(filter(sensemito_fusions, Sequencing_Depth == c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb")), 
       aes(x = Sequencing_Depth, fill = fusionType)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "atypical fusions", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(Library+Cell_Lines~Algorithm) +
  scale_y_log10()+
  labs(fill = "Fusion Type")
