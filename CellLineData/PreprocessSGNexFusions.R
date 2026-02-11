#########################################
# Preprocess Fusion output of SGNex data
#########################################
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))

###########
#LongGF
###########
LongGF_SGNex_Directory<-"/bioinformatics/ryley/Gencode44/SGNex_Fusionv44/LongGF/"
my_files <- list.files(path = LongGF_SGNex_Directory, pattern = "tsv", full.names=TRUE)
all_tsv <-lapply(my_files, read.table)
names(all_tsv) <- gsub(".tsv","",list.files(path = LongGF_SGNex_Directory, pattern = "tsv",full.names = FALSE),fixed = TRUE)
names <- gsub(".tsv","",list.files(path = LongGF_SGNex_Directory, pattern = "tsv",full.names = FALSE),fixed = TRUE)
for(i in names){
  assign(i, read.table(paste0(LongGF_SGNex_Directory, i, ".tsv")))
}
LongGF_SGNex <- do.call(rbind, lapply(my_files, function(filename) {
  read.table(filename) %>%
    mutate(Source = basename(filename))
}))
LongGF_SGNex <- LongGF_SGNex %>% 
  separate(V2, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE) %>% 
  separate(V4, into = c("chromosome1", "breakpoint1"), sep = ":", remove = FALSE) %>% 
  separate(V5, into = c("chromosome2", "breakpoint2"), sep = ":", remove = FALSE) %>%
  filter(chromosome1 %in% standard_chrs)%>%
  filter(chromosome2 %in% standard_chrs)%>%
  filter(V3 >= 2)%>%
  mutate(Algorithm = "LongGF")
  
###################
#Genion
###################
Genion_SGNex_Directory <- "/bioinformatics/ryley/Gencode44/SGNex_Fusionv44/Genion"
myfiles<- list.files(path = Genion_SGNex_Directory, pattern = "_genion$", full.names = TRUE)
Genion_SGNex<- do.call(rbind, lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 0) {
    read.table(filename) %>%
      mutate(Source = sub("_genion$", "", basename(filename)))
  } else {
    NULL
  }
}))
Genion_SGNex<- Genion_SGNex %>% 
  separate(V1, into = c("V1.1", "V1.2", "V1.3"), "::", remove=FALSE) %>% 
  separate(V2, into = c("V2.1", "V2.2", "V2.3"), "::", remove=FALSE) %>% 
  separate(V8, into = c("chr1", "chr2", "chr3"), ";", remove=FALSE) %>%
  separate(chr1, into = c("chr1", "gene1_range"), ":", remove=FALSE)%>%
  separate(chr2, into = c("chr2", "gene2_range"), ":", remove=FALSE)%>%
  separate(chr3, into = c("chr3", "gene3_range"), ":", remove=FALSE)%>%
  filter(V5 >= 2)%>%
  mutate(Algorithm = "Genion")

##############
#JAFFAL
##############
JAF_SGNex_Directory <- "/bioinformatics/ryley/Gencode44/SGNex_Fusionv44/JAFFAL_v43/"
myfiles<-list.files(path = JAF_SGNex_Directory, pattern = "jaffa_results.csv", full.names = TRUE, recursive = TRUE)
JAFFAL_SGNex <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = ",") %>%
    mutate(Source =  sub(".fastq$", "",sample))
}))
JAFFAL_SGNex <- JAFFAL_SGNex %>% 
  separate(fusion.genes, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE) %>%
  filter(chrom1 %in% standard_chrs)%>%
  filter(chrom2 %in% standard_chrs)%>%
  filter(spanning.reads >= 2)%>%
  mutate(Algorithm = "JAFFAL")
  
##############
#JAFFAL_3Gene 
##############
myfiles<-list.files(path = JAF_SGNex_Directory, pattern = ".3gene_summary", full.names = TRUE, recursive = TRUE) 

JAFFAL_SGNex_3Gene <- do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = "\t") %>%
    mutate(Source = sub(".fastq$", "", basename(dirname(filename))))
}))

JAFFAL_SGNex_3Gene <- JAFFAL_SGNex_3Gene %>%
  dplyr::filter(Reads >= 2)%>% 
  separate(Fusion, into = c("Gene1", "Gene2", "Gene3"), sep = ":", remove = FALSE)%>%
  mutate(Algorithm = "JAFFAL")%>%
  mutate(Cell_Lines = case_when(
    grepl("K562", Source, ignore.case = TRUE) ~ "K562",
    grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7"# Keep original value if no match
  ), Sequencing_Depth = case_when(
    grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
    grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
    grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
    grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
    grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
    .default = "Total"), 
  Library = case_when(
    grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
    grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
    .default = "PCR-cDNA"))


##############
#FusionSeeker
##############
FusionSeeker_SGNex_Directory <- "/bioinformatics/ryley/Gencode44/SGNex_Fusionv44/FusionSeeker"

myfiles<-list.files(path = FusionSeeker_SGNex_Directory, pattern = "confident_genefusion.txt", full.names = TRUE, recursive = TRUE)
FusionSeeker_SGNex <- do.call(rbind, lapply(myfiles, function(filename){
  read.table(filename, header = TRUE) %>% mutate(Source = sub("^fusionseeker_", "", basename(dirname(filename))))
}))

FusionSeeker_SGNex <- FusionSeeker_SGNex %>%
  unite("fusionGene", c(Gene1, Gene2), sep="::", remove= FALSE, na.rm = TRUE)%>%
  mutate(across(c(Gene1, Gene2), ~ gsub("\\..*","", .)))%>%
  separate(fusionGene, into= c('Gene1_vID','Gene2_vID'), sep = "::", remove = FALSE) %>%
  unite("fusion.gene.id", c(Gene1, Gene2), sep=":", remove= FALSE, na.rm = TRUE) %>%
  filter(NumSupp >= 2)%>%
  mutate(Algorithm = "FusionSeeker")

#################
# CTAT-LR-Fusion
#################
myfiles<-list.files(path = "/bioinformatics/siyuan/ctat_lr_results/cellline_rasusa", 
                    pattern = "ctat-LR-fusion.fusion_predictions.abridged.tsv", full.names = TRUE, recursive = TRUE)
CTATLR_SGNex<- bind_rows(
  lapply(myfiles, function(filename) {
    if (file.info(filename)$size > 136) {
      read_tsv(filename,
               col_names = TRUE,
               show_col_types = FALSE) %>%
        mutate(Source = basename(dirname(filename)),
               Cell_Lines = case_when(
                 grepl("K562", Source, ignore.case = TRUE) ~ "K562",
                 grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7"
               ), Algorithm = "CTAT-LR-Fusion",
               Sequencing_Depth = case_when(
                 grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
                 grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
                 grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
                 grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
                 grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
                 .default = "Total"),
               Library = case_when(
                 grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
                 grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
                 .default = "PCR-cDNA"))
    } else {
      NULL
    }
  }))
 
CTATLR_SGNex <- CTATLR_SGNex  %>%
  filter(num_LR >= 2)%>% 
  separate(LeftBreakpoint, into = c("chrom1", "base1", "strand1"), ":", remove=FALSE) %>% 
  separate(RightBreakpoint, into = c("chrom2", "base2", "strand2"), ":", remove=FALSE)%>%
  mutate(Algorithm = "CTAT-LR-Fusion")
  
og_CTATLR_SGNex <- CTATLR_SGNex

#################
# GFSeeker
#################
myfiles<-list.files(path = "/bioinformatics/siyuan/gfseeker_results_v2/cellline", 
                    pattern = "final_validated_fusions.csv", full.names = TRUE, recursive = TRUE)
GFSeeker_SGNex<- bind_rows(
  lapply(myfiles, function(filename) {
    if (file.info(filename)$size > 118) {
      read_tsv(filename,
               col_names = TRUE,
               show_col_types = FALSE) %>%
        mutate(Source = basename(dirname(filename)),
               Cell_Lines = case_when(
                 grepl("K562", Source, ignore.case = TRUE) ~ "K562",
                 grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7"
               ), Algorithm = "GFSeeker",
               Sequencing_Depth = case_when(
                 grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
                 grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
                 grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
                 grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
                 grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
                 .default = "Total"),
               Library = case_when(
                 grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
                 grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
                 .default = "PCR-cDNA"))
    } else {
      NULL
    }
  }))

GFSeeker_SGNex <- GFSeeker_SGNex  %>%
  filter("support num" >= 2)%>% 
  rename(gene1_name = "#gene1 name")%>%
  rename(gene2_name = "gene2 name")%>%
  separate("break points", into = c("Breakpoint1", "Breakpoint2"), ";", remove=TRUE) %>% 
  separate(Breakpoint1, into = c("chrom1", "base1"), ":", remove=FALSE)%>%
  mutate(Breakpoint2, gsub("^ " , "", Breakpoint2))%>%
  separate(Breakpoint2, into = c("chrom2", "base2"), ":", remove=FALSE)%>%
  mutate(Algorithm = "GFSeeker")

og_GFSeeker_SGNex <- GFSeeker_SGNex

##############################
# Append Cell lines and Depth
##############################

#Make a column to indicate which cell line the fusion was found in, 
# and column to indicate the sequencing depth and library type

LongGF_SGNex <- LongGF_SGNex %>%
  mutate(Cell_Lines = case_when(
    grepl("K562", Source, ignore.case = TRUE) ~ "K562",
    grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7"
  ),
  Sequencing_Depth = case_when(
    grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
    grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
    grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
    grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
    grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
    .default = "Total"),
  Library = case_when(
    grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
    grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
    .default = "PCR-cDNA"))

Genion_SGNex <- Genion_SGNex%>%
  mutate(Cell_Lines = case_when(
    grepl("K562", Source, ignore.case = TRUE) ~ "K562",
    grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7",# Keep original value if no match
  ),
  Sequencing_Depth = case_when(
    grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
    grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
    grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
    grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
    grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
    .default = "Total"),
  Library = case_when(
    grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
    grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
    .default = "PCR-cDNA"))

FusionSeeker_SGNex <- FusionSeeker_SGNex %>%
  mutate(Cell_Lines = case_when(
    grepl("K562", Source, ignore.case = TRUE) ~ "K562",
    grepl("MCF7", Source, ignore.case = TRUE) ~ "MCF7",# Keep original value if no match
  ),
  Sequencing_Depth = case_when(
    grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
    grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
    grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
    grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
    grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
    .default = "Total"),
  Library = case_when(
    grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
    grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
    .default = "PCR-cDNA"))

JAFFAL_SGNex <- JAFFAL_SGNex%>%
  mutate(Cell_Lines = case_when(
    grepl("K562", sample, ignore.case = TRUE) ~ "K562",
    grepl("MCF7", sample, ignore.case = TRUE) ~ "MCF7"),# Keep original value if no match
    Sequencing_Depth = case_when(
      grepl("1G", Source, ignore.case = TRUE) ~ "1Gb",
      grepl("10G", Source, ignore.case = TRUE) ~ "10Gb",
      grepl("_5G", Source, ignore.case = TRUE) ~ "5Gb",
      grepl("2.5G", Source, ignore.case = TRUE) ~ "2.5Gb",
      grepl("7.5G", Source, ignore.case = TRUE) ~ "7.5Gb",
      .default = "Total"), 
    Library = case_when(
      grepl("RNA", Source, ignore.case = TRUE) ~ "direct-RNA",
      grepl("directcDNA", Source, ignore.case = TRUE) ~ "direct-cDNA",
      .default = "PCR-cDNA"))
 
