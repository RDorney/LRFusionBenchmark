#################################################################
####Huh7 Sequencing Library Fusion Calling Results formatting####
#Author: Ryley Dorney
#Date Feb 2026
#################################################################
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
#Read in data from JAFFA/JAFFAL
################################################################
#Load JAFFAL data into R
JAFFAL_LR_Huh7 <- read.csv("/bioinformatics/ryley/Library_Benchmark/JAFFAL_Huh7_gv43/JAFFAL_Huh7_genCode43/jaffa_results.csv", header = TRUE) %>% 
  mutate(library_type = case_when(grepl('FBA22517|FBA22660|dRNA', sample)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655|dcDNA', sample)  ~ "direct_cDNA",
                                  grepl('FAZ|PCR|FLNCcDNA', sample)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', sample)  ~ "B1", 
                                grepl('FBA22517|FAZ83542|FBA43334|B2', sample) ~ "B2"), 
         Platform = case_when(grepl('Nano', sample)  ~ "ONT",
                              grepl('PB', sample)  ~ "PacBio"),
         Cell_Line = "Huh7") %>%
  mutate(Algorithm = "JAFFAL")
#Load in data from JAFFA direct
Illumina_sequencing_IllHB1<- read.csv("/bioinformatics/ryley/JAFFA_Promethion/Huh7B1/Huh7B1.summary", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",           
         RNA_sample = "B1",
         Platform = "Illumina",
         Algorithm = "JAFFA")
Illumina_sequencing_IllHB2<- read.csv("/bioinformatics/ryley/JAFFA_Promethion/Huh7B2/Huh7B2.summary", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",
         RNA_sample = "B2",
         Platform = "Illumina",
         Algorithm = "JAFFA")
JAFFAL_Illumina_Huh7 <- read.csv("/bioinformatics/ryley/JAFFA_Promethion/jaffa_results.csv", header = TRUE)%>%
  mutate(Algorithm = "JAFFA",
         Platform = "Illumina",
         library_type = "PCR_cDNA") %>%
  mutate(RNA_sample = case_when(grepl('B1', sample)  ~ "B1", 
                                grepl('B2', sample) ~ "B2"))

Huh7_JAFFAL <- rbind(JAFFAL_Illumina_Huh7, JAFFAL_LR_Huh7) %>% #
  separate(fusion.genes, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE) %>%
  filter(chrom1 %in% standard_chrs)%>%
  filter(chrom2 %in% standard_chrs)%>%
  filter(spanning.reads >= 2 |
           spanning.pairs >= 2 |
           (spanning.reads >= 1 & spanning.pairs >= 1)
  )

##load tri gene results
LRTri_gene <- list.files(path="/bioinformatics/ryley/Library_Benchmark/JAFFAL_Huh7_gv43", 
                         pattern = "3gene_summary", full.names = TRUE, recursive = TRUE) 
name <- 0
for (location in LRTri_gene){
  name <- name + 1
  assign(paste0("file_", name), 
         mutate(read.table(location, header = TRUE), 
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

Huh7_JAFFAL_3Gene <- rbind(file_1, file_2, 
                            file_3, file_4, 
                            file_5, file_6,
                            file_7, file_8)%>% 
  separate(Fusion, sep=":" ,into = c("Gene1", "Gene2", "Gene3"), 
           remove = FALSE)

################################################################
#Read in data from Genion
################################################################
Genion_Huh7_Directory <- "/bioinformatics/ryley/Gencode44/Huh7_Library/Genion"
myfiles<- list.files(path = Genion_Huh7_Directory, pattern = "_genion$", full.names = TRUE)
Genion_Huh7<- do.call(rbind, lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 0) {
    read.table(filename) %>%
      mutate(Source = sub("_genion$", "", basename(filename)))
  } else {
    NULL
  }
}))
Genion_Huh7<- Genion_Huh7 %>% 
  separate(V1, into = c("V1.1", "V1.2", "V1.3"), "::", remove=FALSE) %>% 
  separate(V2, into = c("V2.1", "V2.2", "V2.3"), "::", remove=FALSE) %>% 
  separate(V8, into = c("chr1", "chr2", "chr3"), ";", remove=FALSE) %>%
  separate(chr1, into = c("chr1", "gene1_range"), ":", remove=FALSE)%>%
  separate(chr2, into = c("chr2", "gene2_range"), ":", remove=FALSE)%>%
  separate(chr3, into = c("chr3", "gene3_range"), ":", remove=FALSE)%>%
  filter(V5 >= 2)%>%
  mutate(Algorithm = "Genion",
         library_type = case_when(grepl('FBA22517|FBA22660|dRNA', Source)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655|dcDNA', Source)  ~ "direct_cDNA",
                                  grepl('FAZ|PCR|FLNCcDNA', Source)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', Source)  ~ "B1", 
                                grepl('FBA22517|FAZ83542|FBA43334|B2', Source) ~ "B2"), 
         Platform = case_when(grepl('chopper', Source)  ~ "ONT",
                              grepl('PB', Source)  ~ "PacBio"),
         Cell_Line = "Huh7")


################################################################
#Read in data from FusionSeeker
################################################################
FusionSeeker_Huh7_Directory <- "/bioinformatics/ryley/Gencode44/Huh7_Library/FusionSeeker"

myfiles<-list.files(path = FusionSeeker_Huh7_Directory, pattern = "confident_genefusion.txt", full.names = TRUE, recursive = TRUE)
FusionSeeker_Huh7 <- do.call(rbind, lapply(myfiles, function(filename){
  read.table(filename, header = TRUE) %>% 
    mutate(Source = sub("^fusionseeker_", "", basename(dirname(filename))))
}))

FusionSeeker_Huh7 <- FusionSeeker_Huh7 %>%
  unite("fusionGene", c(Gene1, Gene2), sep="::", remove= FALSE, na.rm = TRUE)%>%
  mutate(across(c(Gene1, Gene2), ~ gsub("\\..*","", .)))%>%
  separate(fusionGene, into= c('Gene1_vID','Gene2_vID'), sep = "::", remove = FALSE) %>%
  unite("fusion.gene.id", c(Gene1, Gene2), sep=":", remove= FALSE, na.rm = TRUE) %>%
  filter(NumSupp >= 2)%>%
  mutate(Algorithm = "FusionSeeker",
         library_type = case_when(grepl('FBA22517|FBA22660|dRNA', Source)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655|dcDNA', Source)  ~ "direct_cDNA",
                                  grepl('FAZ|PCR|FLNCcDNA', Source)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', Source)  ~ "B1", 
                                grepl('FBA22517|FAZ83542|FBA43334|B2', Source) ~ "B2"), 
         Platform = case_when(grepl('chopper', Source)  ~ "ONT",
                              grepl('PB', Source)  ~ "PacBio"),
         Cell_Line = "Huh7")

################################################################
#Read in data from LongGF
################################################################
LongGF_Huh7_Directory<-"/bioinformatics/ryley/Gencode44/Huh7_Library/LongGF"
grep_command <- paste0("grep SumGF " , LongGF_Huh7_Directory, "/LongGF.run.*log", " > ", LongGF_Huh7_Directory, "/LongGF_Huh7_summary.tsv")
system(grep_command)

LongGF_Huh7 <- read.table(paste0(LongGF_Huh7_Directory, "/LongGF_Huh7_summary.tsv"))
LongGF_Huh7 <- LongGF_Huh7 %>% 
  mutate(V1 = gsub("/bioinformatics/ryley/Gencode44/Huh7_Library/LongGF/LongGF.run.", "", V1)) %>%
  mutate(V1 =gsub(".log:SumGF", "", V1))%>%
  separate(V2, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE) %>% 
  separate(V4, into = c("chromosome1", "breakpoint1"), sep = ":", remove = FALSE) %>% 
  separate(V5, into = c("chromosome2", "breakpoint2"), sep = ":", remove = FALSE) %>%
  filter(chromosome1 %in% standard_chrs)%>%
  filter(chromosome2 %in% standard_chrs)%>%
  filter(V3 >= 2)%>%
  mutate(Algorithm = "LongGF",
         library_type = case_when(grepl('FBA22517|FBA22660|dRNA', V1)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655|dcDNA', V1)  ~ "direct_cDNA",
                                  grepl('FAZ|PCR|FLNCcDNA', V1)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', V1)  ~ "B1", 
                                grepl('FBA22517|FAZ83542|FBA43334|B2', V1) ~ "B2"), 
         Platform = case_when(grepl('chopper', V1)  ~ "ONT",
                              grepl('PB', V1)  ~ "PacBio"),
         Cell_Line = "Huh7")

################################################################
#Read in data from GFSeeker
################################################################
myfiles<-list.files(path = "/bioinformatics/siyuan/Huh7_results/gfseeker", 
                    pattern = "final_validated_fusions.csv", full.names = TRUE, recursive = TRUE)
GFSeeker_Huh7<- bind_rows(
  lapply(myfiles, function(filename) {
    if (file.info(filename)$size > 118) {
      read_tsv(filename,
               col_names = TRUE,
               show_col_types = FALSE) %>%
        mutate(Source = basename(dirname(filename)), 
               Algorithm = "GFSeeker",
               library_type = case_when(grepl('FBA22517|FBA22660|dRNA', Source)  ~ "direct_RNA",
                                        grepl('FBA43334|FBA62655|dcDNA', Source)  ~ "direct_cDNA",
                                        grepl('FAZ|PCR|FLNCcDNA', Source)  ~ "PCR_cDNA"),
               RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', Source)  ~ "B1", 
                                      grepl('FBA22517|FAZ83542|FBA43334|B2', Source) ~ "B2"), 
               Platform = case_when(grepl('Nano', Source)  ~ "ONT",
                                    grepl('PB', Source)  ~ "PacBio"),
               Cell_Line = "Huh7")
    } else {
      NULL
    }
  }))

GFSeeker_Huh7 <- GFSeeker_Huh7  %>%
  filter("support num" >= 2)%>% 
  rename(gene1_name = "#gene1 name")%>%
  rename(gene2_name = "gene2 name")%>%
  separate("break points", into = c("Breakpoint1", "Breakpoint2"), ";", remove=TRUE) %>% 
  separate(Breakpoint1, into = c("chrom1", "base1"), ":", remove=FALSE)%>%
  separate(Breakpoint2, into = c("chrom2", "base2"), ":", remove=FALSE)

og_GFSeeker_Huh7 <- GFSeeker_Huh7

################################################################
#Read in data from CTAT-LR-Fusion
################################################################
myfiles<-list.files(path = "/bioinformatics/siyuan/Huh7_results/ctat_lr", 
                    pattern = "ctat-LR-fusion.fusion_predictions.abridged.tsv", full.names = TRUE, recursive = TRUE)
CTATLR_Huh7<- bind_rows(
  lapply(myfiles, function(filename) {
    if (file.info(filename)$size > 136) {
      read_tsv(filename,
               col_names = TRUE,
               show_col_types = FALSE) %>%
        mutate(Source = basename(dirname(filename)),
               library_type = case_when(grepl('FBA22517|FBA22660|dRNA', Source)  ~ "direct_RNA",
                                        grepl('FBA43334|FBA62655|dcDNA', Source)  ~ "direct_cDNA",
                                        grepl('FAZ|PCR|FLNCcDNA', Source)  ~ "PCR_cDNA"),
               RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247|B1', Source)  ~ "B1", 
                                      grepl('FBA22517|FAZ83542|FBA43334|B2', Source) ~ "B2"), 
               Platform = case_when(grepl('Nano', Source)  ~ "ONT",
                                    grepl('PB', Source)  ~ "PacBio"),
               Cell_Line = "Huh7")
    } else {
      NULL
    }
  }))

CTATLR_Huh7 <- CTATLR_Huh7  %>%
  filter(num_LR >= 2)%>% 
  separate(LeftBreakpoint, into = c("chrom1", "base1", "strand1"), ":", remove=FALSE) %>% 
  separate(RightBreakpoint, into = c("chrom2", "base2", "strand2"), ":", remove=FALSE)%>%
  mutate(Algorithm = "CTAT-LR-Fusion")

og_CTATLR_Huh7 <- CTATLR_Huh7


