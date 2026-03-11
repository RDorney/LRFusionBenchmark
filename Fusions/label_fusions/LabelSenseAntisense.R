#########################################
# Title: Label Sense Antisense Fusions
# Author: Ryley Dorney
# Date: Mar 2026
#########################################
CTATLR_mitocheck_annotated
Genion_mitocheck_annotated
FusionSeeker_mitocheck_annotated
LongGF_mitocheck_annotated
GFSeeker_mitocheck_annotated
JAFFAL_mitocheck_annotated
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

hits <- findOverlaps(genes, genes, ignore.strand = TRUE)
# Filter hits to find genes on opposite strands that overlap
antisense_lookup <- data.frame(
  gene_a_name = genes$gene_name[queryHits(hits)],
  gene_a_id   = sub("\\..*$", "", genes$gene_id[queryHits(hits)]),
  gene_b_name = genes$gene_name[subjectHits(hits)],
  gene_b_id   = sub("\\..*$", "", genes$gene_id[subjectHits(hits)]),
  strand_a    = as.character(strand(genes)[queryHits(hits)]),
  strand_b    = as.character(strand(genes)[subjectHits(hits)])
) %>%
  filter(strand_a != strand_b) %>% # Only keep antisense
  select(gene_a_name, gene_a_id, gene_b_name, gene_b_id) %>%
  distinct()

##################################
# Alternative function application
#get_antisense_overlaps <- function(gene, genes_gr, id_type = c("gene_name", "gene_id")) {
  
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
########################################
########
#CTATLR
########
CTATLR_sensecheck_annotated <- CTATLR_mitocheck_annotated %>%
  left_join(antisense_lookup, by = c("LeftGene" = "gene_a_name", "RightGene" = "gene_b_name")) %>%
  mutate(
    fusionType = if_else(!is.na(gene_a_id), "Sense-Antisense", fusionType)
  ) %>%
  select(-gene_a_id, -gene_b_id) # Clean up join columns
########
# Genion
########
Genion_sensecheck_annotated <- Genion_mitocheck_annotated %>%
  left_join(
    antisense_lookup, 
    by = c("V1.1" = "gene_a_id", "V1.2" = "gene_b_id")
  ) %>%
  mutate(
    # If a match was found in the lookup table, label it
    fusionType = if_else(!is.na(gene_a_name), "Sense-Antisense", fusionType)
  ) %>%
  # Clean up temporary helper columns
  select(-gene_a_name, -gene_b_name)
##############
# FusionSeeker
##############
FusionSeeker_sensecheck_annotated <- FusionSeeker_mitocheck_annotated %>%
  left_join(
    antisense_lookup, 
    by = c("Gene1" = "gene_a_id", "Gene2" = "gene_b_id")
  ) %>%
  mutate(
    # If a match was found in the lookup table, label it
    fusionType = if_else(!is.na(gene_a_name), "Sense-Antisense", fusionType)
  ) %>%
  # Clean up temporary helper columns
  select(-gene_a_name, -gene_b_name)

#########
# LongGF
#########
LongGF_sensecheck_annotated <- LongGF_mitocheck_annotated %>%
  left_join(antisense_lookup, by = c("Gene1" = "gene_a_name", "Gene2" = "gene_b_name")) %>%
  mutate(
    fusionType = if_else(!is.na(gene_a_id), "Sense-Antisense", fusionType)
  ) %>%
  select(-gene_a_id, -gene_b_id) # Clean up join columns
##########
# GFSeeker
##########
GFSeeker_sensecheck_annotated <- GFSeeker_mitocheck_annotated %>%
  left_join(antisense_lookup, by = c("gene1_name" = "gene_a_name", "gene2_name" = "gene_b_name")) %>%
  mutate(
    fusionType = if_else(!is.na(gene_a_id), "Sense-Antisense", fusionType)
  ) %>%
  select(-gene_a_id, -gene_b_id) # Clean up join columns

########################
# JAFFAL/JAFFA-Direct
########################
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

hits <- findOverlaps(genes, genes, ignore.strand = TRUE)
# Filter hits to find genes on opposite strands that overlap
antisense_lookup <- data.frame(
  gene_a_name = genes$gene_name[queryHits(hits)],
  gene_a_id   = sub("\\..*$", "", genes$gene_id[queryHits(hits)]),
  gene_b_name = genes$gene_name[subjectHits(hits)],
  gene_b_id   = sub("\\..*$", "", genes$gene_id[subjectHits(hits)]),
  strand_a    = as.character(strand(genes)[queryHits(hits)]),
  strand_b    = as.character(strand(genes)[subjectHits(hits)])
) %>%
  filter(strand_a != strand_b) %>% # Only keep antisense
  select(gene_a_name, gene_a_id, gene_b_name, gene_b_id) %>%
  distinct()


JAFFAL_sensecheck_annotated <- JAFFAL_mitocheck_annotated %>%
  left_join(antisense_lookup, by = c("Gene1" = "gene_a_name", "Gene2" = "gene_b_name")) %>%
  mutate(
    fusionType = if_else(!is.na(gene_a_id), "Sense-Antisense", fusionType)
  ) %>%
  select(-gene_a_id, -gene_b_id) # Clean up join columns

#############################
#Plot Sense-Antisense Fusions
#############################
sensecheck_fusions <- rbind(dplyr::select(CTATLR_sensecheck_annotated, 
                                         c("RNA_sample", "Platform", "Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")),
                           dplyr::select(Genion_sensecheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")),
                           dplyr::select(LongGF_sensecheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")), 
                           dplyr::select(FusionSeeker_sensecheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")),
                           dplyr::select(GFSeeker_sensecheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")),
                           dplyr::select(JAFFAL_sensecheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType", "full_label")))

sensecheck_fusions$Platform <- factor(sensecheck_fusions$Platform, levels = c("Illumina", "PacBio", "ONT")) 

sensecheck_fusions$library_type <- factor(sensecheck_fusions$library_type, levels = c("PCR_cDNA", "direct_cDNA", "direct_RNA")) 

sensecheck_fusions$full_label <- factor(sensecheck_fusions$full_label, levels = c("Illumina.PCR_cDNA.B1", "Illumina.PCR_cDNA.B2", 
                                                                                "PacBio.PCR_cDNA.B1", "PacBio.PCR_cDNA.B2", 
                                                                                "ONT.PCR_cDNA.B1","ONT.PCR_cDNA.B2", 
                                                                                "ONT.direct_cDNA.B1", "ONT.direct_cDNA.B2", 
                                                                                "ONT.direct_RNA.B1", "ONT.direct_RNA.B2"))


ggplot(dplyr::filter(sensecheck_fusions, fusionType== "Sense-Antisense"), 
       aes(x = RNA_sample, colour = Algorithm)) +
  geom_point(stat = "count") +
  theme_minimal() +
  labs(title = "Sense-Antisense", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  labs(fill = "Algorithm")+
  scale_y_log10()+
  scale_colour_manual(values = Alg_colour_map)

ggplot(dplyr::filter(sensecheck_fusions, fusionType== "Sense-Antisense"), 
       aes(x = RNA_sample, colour = full_label, shape = Algorithm)) +
  geom_point(stat = "count") +
  theme_minimal() +
  labs(title = "Sense-Antisense", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  #labs(fill = "Fusion Type")+
  scale_y_log10()+
  scale_colour_manual(values = platformlibsamp_colourmap)
