####Huh7 Sequencing Library JAFFAL results formatting####
#Author: Ryley Dorney
#Date Dec 30-25
#####
#Load R libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)
library(ComplexUpset)
library(biomaRt)
#Load JAFFAL data into R
LR_sequencing <- read.csv("/bioinformatics/ryley/Library_Benchmark/JAFFAL_Huh7_gv43/JAFFAL_Huh7_genCode43/jaffa_results.csv", header = TRUE) %>% 
  mutate(library_type = case_when(grepl('FBA22517|FBA22660|dRNA', sample)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655|dcDNA', sample)  ~ "direct_cDNA",
                                  grepl('FAZ|PCR|FLNCcDNA', sample)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', sample)  ~ "B1", 
                                grepl('FBA22517|FAZ83542|FBA43334|B2', sample) ~ "B2"), 
         Platform = case_when(grepl('Nano', sample)  ~ "ONT",
                              grepl('PB', sample)  ~ "PacBio"),
         Cell_Line = "Huh7")

Illumina_sequencing_IllHB1<- read.csv("/bioinformatics/ryley/Library_Benchmark/JAFFA_direct_gencode43/Huh7B1/jaffa_results_Huh7B1.csv", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",           
         RNA_sample = "B1",
         Platform = "Illumina")
Illumina_sequencing_IllHB2<- read.csv("/bioinformatics/ryley/Library_Benchmark/JAFFA_direct_gencode43/Huh7B2/jaffa_results_Huh7B2.csv", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",
         RNA_sample = "B2",
         Platform = "Illumina")
Illumina_sequencing <- rbind(Illumina_sequencing_IllHB1, Illumina_sequencing_IllHB2)


Huh7_JAFFAL <- rbind(Illumina_sequencing, LR_sequencing)

##load tri gene results
LRTri_gene <- list.files(path="/bioinformatics/ryley/Library_Benchmark/JAFFAL_Huh7_gv43", pattern = "3gene_summary", full.names = TRUE, recursive = TRUE) 
name <- 0
for (location in LRTri_gene){
  name <- name + 1
  assign(paste0("file_", name), mutate(read.table(location, header = TRUE), 
                                       sample = location, 
                                       library_type = case_when(grepl('FBA22517|FBA22660|dRNA', sample)  ~ "direct_RNA",
                                                                grepl('FBA43334|FBA62655|dcDNA', sample)  ~ "direct_cDNA",
                                                                grepl('FAZ|PCR|FLNCcDNA', sample)  ~ "PCR_cDNA"),
                                       RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', sample)  ~ "B1", 
                                                              grepl('FBA22517|FAZ83542|FBA43334|B2', sample) ~ "B2"), 
                                       Platform = case_when(grepl('Nano', sample)  ~ "ONT",
                                                            grepl('PB', sample)  ~ "PacBio"),
                                       Cell_Line = "Huh7"))
} 

LRsequencing_3Gene <- rbind(file_1, file_2, 
                            file_3, file_4, 
                            file_5, file_6,
                            file_7, file_8)
########################################################
#filtering out rRNA fusions and low read support fusions
obvious_library_chimeras <- Huh7_JAFFAL %>%
  filter(!grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom1) | #non-genomic chromosomes and scaffolds
           !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom2)) 
lowread_Huh7_JAFFAL <- Huh7_JAFFAL %>% filter((spanning.pairs==0 & spanning.reads==1) |(spanning.pairs==1 & spanning.reads==0)) #fusions supported by only one spanning read or one spanning pair


filtered_Huh7_JAFFAL <- Huh7_JAFFAL %>%
  filter(grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom1) &
           grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom2)) %>%
  filter((spanning.pairs>=0 & spanning.reads>=2) |(spanning.pairs>=2 & spanning.reads>=0)) %>% 
  mutate(fusiontype = case_when(chrom1 == chrom2 ~"Intra-chromosomal", chrom1 != chrom2 ~ "Inter-chromosomal"))

filtered_Huh7_JAFFAL_3Gene <- Huh7_JAFFAL_3Gene %>%
  filter(Reads>=2)
##################################
#Re-annotate fusions with geneIDs#
##################################
#Load Ensembl
ensembl_current <-  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
#the below ensembl versions best correlate to gencode version 43
ensemblv109 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 109)
ensemblv110 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 110)
ensembljul2023 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://jul2023.archive.ensembl.org")

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
#Create library type, platform and RNA sample as different factors
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


# Define colors for other platforms
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