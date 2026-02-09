####Huh7 Sequencing Library Fusion Calling Results formatting####
#Author: Ryley Dorney
#Date Dec 30-25
#####
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
################################################################
#Load R libraries
################################################################
library(dplyr)
library(tidyr)

library(biomaRt)
library(AnnotationHub)
library(HGNChelper)

################################################################
#Re-annotate fusions with geneIDs#
################################################################
#Load Ensembl the version that bests correlate to gencode version 44
ensemblv110 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 110)

#the current ensembl version incase any GeneIDs have since been assigned names or vice versa
ensembl_current <-  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

sep_fil_Huh7_JAF <- filtered_Huh7_JAFFAL %>% separate(fusion.genes, sep=":" ,into = c("gene1", "gene2"), remove = FALSE)
filtered_Huh7_JAFFAL_3Gene  %>% separate(fusion.genes, sep=":" ,into = c("gene1", "gene2", "gene3"), remove = FALSE)

GeneIDs<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "external_gene_name",
                        values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                        mart = ensembljul2023),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "external_gene_name",
                        values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                        mart = ensemblv109),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "external_gene_name",
                        values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                        mart = ensemblv110),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                     filters = "external_gene_name",
                     values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                     mart = ensembl_current)) %>% unique()
GeneNames<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                       mart = ensembljul2023),
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                       mart = ensemblv109),
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                       mart = ensemblv110),
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = c(sep_fil_Huh7_JAF$gene1, sep_fil_Huh7_JAF$gene2),
                       mart = ensembl_current)) %>% unique()
EnsemblInfo <- rbind(GeneIDs, GeneNames)
################################################################
#Create library type, platform and RNA sample as different factors
################################################################
Huh7_JAFFAL$library_type <- factor(
  Huh7_JAFFAL$library_type,
  levels = c("direct_RNA","direct_cDNA", "PCR_cDNA")
)
Huh7_JAFFAL$Platform <- factor(
  Huh7_JAFFAL$Platform,
  levels = c("ONT", "PacBio", "Illumina")
)
Huh7_JAFFAL$RNA_sample <- factor(
  Huh7_JAFFAL$RNA_sample,
  levels = c("B2", "B1"))

################################################################
# Define colors for other platforms
################################################################
platform_colours <- c(
  "Illumina.PCR_cDNA.B1" = "#FFA500",
  "Illumina.PCR_cDNA.B2" = "#CC8400",
  "PacBio.PCR_cDNA.B1"   = "#C80972",
  "PacBio.PCR_cDNA.B2"   = "#8B064F",
  "ONT.direct_cDNA.B1" = "#009ACD",
  "ONT.direct_cDNA.B2" = "#00688B",
  "ONT.direct_RNA.B1"  = "#00CDCD",
  "ONT.direct_RNA.B2"  = "#008B8B",
  "ONT.PCR_cDNA.B1" = "#1f77b4",
  "ONT.PCR_cDNA.B2" = "#15496e"
)
# Merge colors
color_map <- c(platform_colours)