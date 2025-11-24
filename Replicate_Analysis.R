#Replicate Analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
library(biomaRt)

#Load biomaRt####
ensembl_current <-  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
ensemblv80 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 80)
ensemblv109 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 109)
ensemblv110 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 110)
ensembljul2023 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://jul2023.archive.ensembl.org")
#read simulated fusions####
simulated_fusions <- read.table(file = "~/benchmark_data_analysis/Jul23_benchmark_fusions.txt", header=TRUE, sep= "\t")
write_tsv(simulated_fusions, file="~/simulated_fusion_info.tsv")

Simulated_Fusion_Info <- simulated_fusions %>% separate(fusionGene, into= c('V1','V2','V3'), sep = "-", remove = FALSE) %>% mutate(across(fusionGene, ~ gsub("-", "::", .))) %>% mutate(geneName_ID = geneName) %>% mutate(across(c(V1, V2, V3, geneName), ~ gsub("\\..*","", .)))%>% separate(fusionGene, into= c('Gene1_vID','Gene2_vID','Gene3_vID'), sep = "::", remove = FALSE)

Gene1_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "ensembl_gene_id",
                        values = (Simulated_Fusion_Info$geneName),
                        mart = ensembljul2023),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "ensembl_gene_id",
                        values = (Simulated_Fusion_Info$geneName),
                        mart = ensemblv109),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "ensembl_gene_id",
                        values = (Simulated_Fusion_Info$geneName),
                        mart = ensemblv110)) %>% unique()
External_Gene_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = "external_gene_name",
                                values = unique(Gene1_Name$external_gene_name),
                                mart = ensembljul2023), 
                          getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = "external_gene_name",
                                values = unique(Gene1_Name$external_gene_name),
                                mart = ensemblv109), 
                          getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = "external_gene_name",
                                values = unique(Gene1_Name$external_gene_name),
                                mart = ensemblv110)) %>% unique()
External_Gene_Name$chromosome_name <- paste0("chr",External_Gene_Name$chromosome_name)

Alternative_Names_Gene_Fusion <- rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                                             filters = "external_gene_name",
                                             values = (External_Gene_Name$external_gene_name),
                                             mart = ensembljul2023),
                                       getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                                             filters = "external_gene_name",
                                             values = (External_Gene_Name$external_gene_name),
                                             mart = ensemblv109),
                                       getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                                             filters = "external_gene_name",
                                             values = (External_Gene_Name$external_gene_name),
                                             mart = ensemblv110)) %>% unique()
Alternative_Names_Gene_Fusion$chromosome_name <- paste0("chr", Alternative_Names_Gene_Fusion$chromosome_name)

Sim_Info_Names <- left_join(Simulated_Fusion_Info, Alternative_Names_Gene_Fusion[c(1,2)], by = c("geneName"="ensembl_gene_id"))%>% 
  left_join(Alternative_Names_Gene_Fusion, by = c("external_gene_name", "chrom" = "chromosome_name")) %>% mutate(
    alternative_ID = coalesce(ensembl_gene_id, geneName)) %>% select(-c(ensembl_gene_id))

Simulated_Fusion_Info_2<-Sim_Info_Names %>%
  mutate(Gene1_alternative_ID = map(Gene1_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene2_alternative_ID = map(Gene2_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene3_alternative_ID = map(Gene3_vID, ~ if (!is.na(.x)) alternative_ID[which(geneName_ID == .x)] else NA_character_)) %>%
  unnest_longer(Gene1_alternative_ID) %>%
  unnest_longer(Gene2_alternative_ID) %>%
  unnest_longer(Gene3_alternative_ID, keep_empty = TRUE)%>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID, Gene3_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>%
  select(-c(alternative_ID)) %>% unique() %>%
  unite("original.fusion.gene.id", c(V1, V2), sep=":", remove= FALSE, na.rm = TRUE)

length(unique(Simulated_Fusion_Info_2$fusionGene))

Simulated_Fusion_Info_3<-Sim_Info_Names %>%
  mutate(Gene1_alternative_ID = map(Gene1_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene2_alternative_ID = map(Gene2_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene3_alternative_ID = map(Gene3_vID, ~ if (!is.na(.x)) alternative_ID[which(geneName_ID == .x)] else NA_character_)) %>%
  unnest_longer(Gene1_alternative_ID) %>%
  unnest_longer(Gene2_alternative_ID) %>%
  unnest_longer(Gene3_alternative_ID, keep_empty = TRUE)%>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID, Gene3_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>%
  unique() %>%
  unite("original.fusion.gene.id", c(V1, V2), sep=":", remove= FALSE, na.rm = TRUE)

#Read in data from LongGF####
system("grep SumGF /bioinformatics/ryley/Algorithm_Benchmark/Replicates/LongGF/LongGF*Spiked*.log > /home/ryleyd/LongReadFusionCallerBenchmark/Replicates/LongGF_Spiked.tsv")

#Process LongGF data
LongGF_Replicate_Sim <-  read.table("/home/ryleyd//Replicates/LongGF_Spiked.tsv")%>%
  mutate(Source = basename(V1)) %>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))%>% 
  separate(V2, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE) %>% 
  separate(V4, into = c("chromosome1", "breakpoint1"), sep = ":", remove = FALSE) %>% 
  separate(V5, into = c("chromosome2", "breakpoint2"), sep = ":", remove = FALSE) #V2 = fusion gene name :
ogLongGF_Replicate <- LongGF_Replicate_Sim

Gene_Name<-rbind (getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                        filters = "external_gene_name",
                        values = unique(c(LongGF_Replicate_Sim$Gene1, LongGF_Replicate_Sim$Gene2)),
                        mart = ensembljul2023),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                        filters = "external_gene_name",
                        values = unique(c(LongGF_Replicate_Sim$Gene1, LongGF_Replicate_Sim$Gene2)),
                        mart = ensemblv109)) %>% unique()
Gene_Name$chromosome_name <- paste0("chr", Gene_Name$chromosome_name)

LongGF_Replicate_Sim_ensembl <- LongGF_Replicate_Sim %>%
  left_join(Gene_Name, by = c('Gene1' = 'external_gene_name')) %>%
  left_join(Gene_Name, by = c('Gene2' = 'external_gene_name')) %>%
  mutate(
    ensembl_gene1 = ifelse(is.na(ensembl_gene_id.x) | nchar(ensembl_gene_id.x) == 0, Gene1, ensembl_gene_id.x),
    ensembl_gene2 = ifelse(is.na(ensembl_gene_id.y) | nchar(ensembl_gene_id.y) == 0, Gene2, ensembl_gene_id.y)
  ) %>%
  unite("fusion.gene.id", c(ensembl_gene1, ensembl_gene2), sep = ":", remove = FALSE, na.rm = TRUE) %>%
  # Create a filtering condition using case_when
  filter(case_when(
    is.na(chromosome_name.x) ~ TRUE, # Keep if chromosome_name.x is NA
    chromosome_name.x == chromosome1 ~ TRUE, # Keep if chromosome names match
    chromosome_name.x == "chrMT" & chromosome1 == "chrM" ~ TRUE, # Allow chrMT/chrM equivalence
    grepl("PATCH|HSCHR", chromosome_name.x) & !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$", chromosome1) ~ TRUE, 
    TRUE ~ FALSE # Remove rows that don't meet the above conditions
  )) %>%
  filter(case_when(
    is.na(chromosome_name.y) ~ TRUE,
    chromosome_name.y == chromosome2 ~ TRUE,
    chromosome_name.y == "chrMT" & chromosome2 == "chrM" ~ TRUE,
    grepl("PATCH|HSCHR", chromosome_name.y) & !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$", chromosome2) ~ TRUE,
    TRUE ~ FALSE # Remove rows that don't meet the above conditions
  )) 

Annot_LongGF_Sim_Replicate <- LongGF_Replicate_Sim_ensembl %>% left_join(Simulated_Fusion_Info_2, by = c("fusion.gene.id"), relationship = "many-to-many")  
Annot_LongGF_Sim_Replicate <- Annot_LongGF_Sim_Replicate %>%
  group_by(V2.x) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA")) 

Annot_LongGF_Sim_Replicate <- LongGF_Replicate_Sim_ensembl %>% left_join(Simulated_Fusion_Info_2, by = c("fusion.gene.id"), relationship = "many-to-many") %>%
  group_by(V2.x) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA"))  
Annot_LongGF_Sim_Replicate$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for truncated_tri_fusions
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
            str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
      return("truncated_tri_fusion")
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1))|
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
      return("reverse_order:truncated_tri_fusion") 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
      return(paste("reverse_order:", matching_row[1])) 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
    } else if ((grepl("chrM:", chr1, ignore.case = TRUE) & !grepl("chrM:", chr2, ignore.case = TRUE)) | (!grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE))){
      return("false_fusion:mitochondrial_genomic") 
    }else if (grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE)){
      return("false_fusion:mitochondrial") 
    }else if ((gene_name1 == gene_name2) & (chr1 != chr2)){
      return("false_fusion:self_misalignment") 
    }else {
      return("false_fusion")
    }
  }
  return(current_type)
}, Annot_LongGF_Sim_Replicate$ensembl_gene1, Annot_LongGF_Sim_Replicate$ensembl_gene2, Annot_LongGF_Sim_Replicate$fusionType, Annot_LongGF_Sim_Replicate$chromosome1, Annot_LongGF_Sim_Replicate$chromosome2, Annot_LongGF_Sim_Replicate$Gene1, Annot_LongGF_Sim_Replicate$Gene2)

#FusionSeeker Replicate Analysis####
myfiles<-list.files(path = "/bioinformatics/ryley/Algorithm_Benchmark/Replicates/FusionSeeker", pattern = "confident_genefusion.txt", full.names = TRUE, recursive = TRUE)
FusionSeeker_Sim_Replicate <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE) %>%
    mutate(Source = basename(dirname(filename)))%>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
}))
ogFusionSeeker_Replicate <- do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE) %>%
    mutate(Source = basename(dirname(filename)))%>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
}))

FusionSeeker_Sim_Replicate <- FusionSeeker_Sim_Replicate %>%
  unite("fusionGene", c(Gene1, Gene2), sep="::", remove= FALSE, na.rm = TRUE)%>%
  mutate(across(c(Gene1, Gene2), ~ gsub("\\..*","", .)))%>%
  separate(fusionGene, into= c('Gene1_vID','Gene2_vID'), sep = "::", remove = FALSE) %>%
  unite("original.fusion.gene.id", c(Gene1, Gene2), sep=":", remove= FALSE, na.rm = TRUE)
PAR_Y_problem1 <- FusionSeeker_Sim_Replicate %>%
  filter(grepl("PAR_Y", Gene1_vID, ignore.case = TRUE)) %>% select(Gene1) 
PAR_Y_problem2 <-  FusionSeeker_Sim_Replicate %>%
  filter(grepl("PAR_Y", Gene2_vID, ignore.case = TRUE)) %>% select(Gene2)

Gene_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = (unique(c(PAR_Y_problem1$Gene1, PAR_Y_problem2$Gene2))),
                       mart = ensembljul2023),
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = (unique(c(PAR_Y_problem1$Gene1, PAR_Y_problem2$Gene2))),
                       mart = ensemblv109),
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = (unique(c(PAR_Y_problem1$Gene1, PAR_Y_problem2$Gene2))),
                       mart = ensemblv110)) %>% unique()
External_Gene_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = c("chromosome_name", "external_gene_name"),
                                values = list("Y", unique(Gene_Name$external_gene_name)),
                                mart = ensembljul2023), 
                          getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = c("chromosome_name", "external_gene_name"),
                                values = list("Y", unique(Gene_Name$external_gene_name)),
                                mart = ensemblv109), 
                          getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version",  "chromosome_name"),
                                filters = c("chromosome_name", "external_gene_name"),
                                values = list("Y", unique(Gene_Name$external_gene_name)),
                                mart = ensemblv110)) %>% unique()
External_Gene_Name$chromosome_name <- paste0("chr", External_Gene_Name$chromosome_name)
FusionSeeker_Sim_Replicate_ensembl <- FusionSeeker_Sim_Replicate %>% left_join(select((right_join(Gene_Name, External_Gene_Name, by= "external_gene_name")), c(ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id_version, chromosome_name)), by= c('Gene1'='ensembl_gene_id.x', 'Chrom1'='chromosome_name')) %>% 
  mutate(Gene1_alternative_ID= coalesce(ensembl_gene_id.y, Gene1), Gene1_alternative_vID= coalesce(ensembl_gene_id_version, Gene1_vID)) %>%
  select(-c(ensembl_gene_id.y, ensembl_gene_id_version))%>% 
  left_join(select((right_join(Gene_Name, External_Gene_Name, by= "external_gene_name")), c(ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id_version, chromosome_name)), by= c('Gene2'='ensembl_gene_id.x', 'Chrom2'='chromosome_name')) %>% 
  mutate(Gene2_alternative_ID= coalesce(ensembl_gene_id.y, Gene2), Gene2_alternative_vID= coalesce(ensembl_gene_id_version, Gene2_vID)) %>%
  select(-c(ensembl_gene_id.y, ensembl_gene_id_version)) %>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>% unique() 

setdiff(unique(FusionSeeker_Sim_Replicate_ensembl$fusion.genes), unique(FusionSeeker_Sim_Replicate$fusion.genes))

Annot_FusionSeeker_Sim_Replicate <- FusionSeeker_Sim_Replicate_ensembl%>%
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id')%>%
  group_by(fusionGene.x) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA"))

Annot_FusionSeeker_Sim_Replicate$fusionType <- mapply(function(g1, g2, Gen1, Gen2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
            str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
      return("truncated_tri_fusion")
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)) |
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
      return("reverse_order:truncated_tri_fusion") 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
      return(paste("reverse_order:", matching_row[1])) 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste0("chromosomal_misalignment:", matching_row[1])) 
    } else if ((grepl("chrM", chr1, ignore.case = TRUE) & !grepl("chrM", chr2, ignore.case = TRUE)) | (!grepl("chrM", chr1, ignore.case = TRUE) & grepl("chrM", chr2, ignore.case = TRUE))){
      return("false_fusion:mitochondrial_genomic") 
    }else if (grepl("chrM", chr1, ignore.case = TRUE) & grepl("chrM", chr2, ignore.case = TRUE)){
      return("false_fusion:mitochondrial") 
    }else if (((g1 == g2)|(Gen1 == Gen2)) & (chr1 != chr2)){
      return("false_fusion:self_misalignment") 
    }else {
      return("false_fusion")
    }
  }
  return(current_type)
}, Annot_FusionSeeker_Sim_Replicate$Gene1_alternative_ID.x, Annot_FusionSeeker_Sim_Replicate$Gene2_alternative_ID.x,
Annot_FusionSeeker_Sim_Replicate$Gene1, Annot_FusionSeeker_Sim_Replicate$Gene2, 
Annot_FusionSeeker_Sim_Replicate$fusionType, Annot_FusionSeeker_Sim_Replicate$Chrom1, Annot_FusionSeeker_Sim_Replicate$Chrom2)
Annot_FusionSeeker_Sim <- Annot_FusionSeeker_Sim %>%
  mutate(depth = case_when(
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ), Replicate_Seed = case_when(
    grepl("replicate_1.", sample) ~ "1",
    grepl("replicate_2", sample) ~ "2",
    grepl("replicate_3", sample) ~ "3",
    grepl("replicate_4", sample) ~ "4",
    grepl("replicate_5", sample) ~ "5",
    grepl("replicate_6", sample) ~ "6",
    grepl("replicate_7", sample) ~ "7",
    grepl("replicate_8", sample) ~ "8",
    grepl("replicate_9", sample) ~ "9",
    grepl("replicate_10", sample) ~ "10"
  ))

#JAFFAL####
#read in data
JAFFAL_Sim_Replicate <- read.table("/bioinformatics/ryley/Algorithm_Benchmark/Replicates/jaffal_results_replicates/jaffa_results.csv", header = TRUE, sep = ",") %>%
     mutate(control = ifelse(grepl("Spiked", sample), "positive", "negative"))%>% 
  separate(fusion.genes, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE)
ogJAFFAL_replicate <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = ",") %>%
    mutate(Source = basename(dirname(filename))) %>% mutate(control = ifelse(grepl("Spiked", sample), "positive", "negative"))
}))

# fusion.genes = fusion gene name : 
Gene_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = unique(c(JAFFAL_Sim_Replicate$Gene1, JAFFAL_Sim_Replicate$Gene2)),
                       mart = ensembljul2023) , getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                                                      filters = "external_gene_name",
                                                      values = unique(c(JAFFAL_Sim_Replicate$Gene1, JAFFAL_Sim_Replicate$Gene2)),
                                                      mart = ensemblv109)) %>% unique()
Gene_Name$chromosome_name <- paste0("chr", Gene_Name$chromosome_name)

JAFFAL_Sim_Replicate_ensembl <- left_join(JAFFAL_Sim_Replicate, Gene_Name, by= c('Gene1'='external_gene_name')) %>% 
  left_join(Gene_Name, by= c('Gene2'='external_gene_name'))  %>%
  mutate(ensembl_gene1 = ifelse(is.na(ensembl_gene_id.x) |nchar(ensembl_gene_id.x) == 0, Gene1, ensembl_gene_id.x), 
         ensembl_gene2 = ifelse(is.na(ensembl_gene_id.y) |nchar(ensembl_gene_id.y) == 0, Gene2, ensembl_gene_id.y)) %>%
  unite("fusion.gene.id", c(ensembl_gene1, ensembl_gene2), sep=":", remove= FALSE, na.rm = TRUE) %>%
  # Create a filtering condition using case_when
  filter(case_when(
    is.na(chromosome_name.x) ~ TRUE, # Keep if chromosome_name.x is NA
    chromosome_name.x == chrom1 ~ TRUE, # Keep if chromosome names match
    chromosome_name.x == "chrMT" & chrom1 == "chrM" ~ TRUE, # Allow chrMT/chrM equivalence
    grepl("PATCH|HSCHR", chromosome_name.x) & !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$", chrom1) ~ TRUE, 
    TRUE ~ FALSE # Remove rows that don't meet the above conditions
  )) %>%
  filter(case_when(
    is.na(chromosome_name.y) ~ TRUE,
    chromosome_name.y == chrom2 ~ TRUE,
    chromosome_name.y == "chrMT" & chrom2 == "chrM" ~ TRUE,
    grepl("PATCH|HSCHR", chromosome_name.y) & !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$", chrom2) ~ TRUE,
    TRUE ~ FALSE # Remove rows that don't meet the above conditions
  )) %>%  
  select(-c("ensembl_gene_id.x", "chromosome_name.x", "ensembl_gene_id.y", "chromosome_name.y")) %>% unique()

setdiff(unique(JAFFAL_Sim_Replicate_ensembl$fusion.genes), unique(ogJAFFAL$fusion.genes))

#Annotate as known or false positive
Annot_JAFFAL_Sim_Replicate <- JAFFAL_Sim_Replicate_ensembl %>% 
  left_join(Simulated_Fusion_Info_2, by = c('fusion.gene.id')) %>% unique()%>%
  group_by(fusion.genes) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA")) 
# Check if ensembl_gene1 and ensembl_gene2 are contained in fusion.gene.id
Annot_JAFFAL_Sim_Replicate$fusionType <- mapply(function(g1, g2, Gen1, Gen2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2))|
            str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
      return("truncated_tri_fusion")
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
      return(paste("reverse_order:", matching_row[1])) 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1))|
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
      return("reverse_order:truncated_tri_fusion") 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
    }else if ((grepl("chrM:", chr1, ignore.case = TRUE) & !grepl("chrM:", chr2, ignore.case = TRUE)) | (!grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE))){
      return("false_fusion:mitochondrial_genomic") 
    }else if (grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE)){
      return("false_fusion:mitochondrial") 
    }else if (((g1 == g2)|(Gen1 == Gen2)) & (chr1 != chr2)){
      return("false_fusion:self_misalignment") 
    }else {
      return("false_fusion")
    }
  }
  return(current_type)
}, Annot_JAFFAL_Sim_Replicate$ensembl_gene1, Annot_JAFFAL_Sim_Replicate$ensembl_gene2, Annot_JAFFAL_Sim_Replicate$Gene1, Annot_JAFFAL_Sim_Replicate$Gene2, Annot_JAFFAL_Sim_Replicate$fusionType, Annot_JAFFAL_Sim_Replicate$chrom1, Annot_JAFFAL_Sim_Replicate$chrom2)

Annot_JAFFAL_Sim_Replicate <- Annot_JAFFAL_Sim_Replicate%>%
  mutate(depth = case_when(
    grepl("10GB", sample) ~ "10GB",
    grepl("1GB", sample) ~ "1GB"# Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", sample) ~ "95%",
    grepl("90", sample) ~ "90%",
    grepl("85", sample) ~ "85%" # Keep original value if no match
  ), Replicate_Seed = case_when(
    grepl("replicate_1.", sample) ~ "1",
    grepl("replicate_2", sample) ~ "2",
    grepl("replicate_3", sample) ~ "3",
    grepl("replicate_4", sample) ~ "4",
    grepl("replicate_5", sample) ~ "5",
    grepl("replicate_6", sample) ~ "6",
    grepl("replicate_7", sample) ~ "7",
    grepl("replicate_8", sample) ~ "8",
    grepl("replicate_9", sample) ~ "9",
    grepl("replicate_10", sample) ~ "10"
  ))
#JAFFAL 3 Gene ####
myfiles<-list.files(path = "/bioinformatics/ryley/Algorithm_Benchmark/Replicates/jaffal_results_replicates/", pattern = ".3gene_summary", full.names = TRUE, recursive = TRUE)
JAFFAL_3Gene_Sim_Replicate <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = "\t") %>%
    mutate(Source = basename(dirname(filename))) %>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
})) %>% separate(Fusion, into = c("Gene1", "Gene2", "Gene3"), sep = ":", remove = FALSE)
Gene_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = unique(c(JAFFAL_3Gene_Sim_Replicate$Gene1, JAFFAL_3Gene_Sim_Replicate$Gene2, JAFFAL_3Gene_Sim_Replicate$Gene3)),
                       mart = ensembljul2023) ,
                 getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = unique(c(JAFFAL_3Gene_Sim_Replicate$Gene1, JAFFAL_3Gene_Sim_Replicate$Gene2, JAFFAL_3Gene_Sim_Replicate$Gene3)),
                       mart = ensemblv109))

JAFFAL_3Gene_Sim_Replicate <- left_join(JAFFAL_3Gene_Sim_Replicate, Gene_Name, by= c('Gene1'='external_gene_name')) %>%
  left_join(Gene_Name, by= c('Gene2'='external_gene_name')) %>%
  left_join(Gene_Name, by= c('Gene3'='external_gene_name')) %>% unique() %>%
  mutate(ensembl_gene1 = ifelse(is.na(ensembl_gene_id.x) |nchar(ensembl_gene_id.x) == 0, Gene1, ensembl_gene_id.x), 
         ensembl_gene2 = ifelse(is.na(ensembl_gene_id.y) |nchar(ensembl_gene_id.y) == 0, Gene2, ensembl_gene_id.y), 
         ensembl_gene3 = ifelse(is.na(ensembl_gene_id) |nchar(ensembl_gene_id) == 0, Gene3, ensembl_gene_id)) %>%
  unite("fusion.gene.id", c(ensembl_gene1, ensembl_gene2, ensembl_gene3), sep=":", remove= FALSE, na.rm = TRUE)

Annot_JAFFAL_3Gene_Sim_Replicate <- JAFFAL_3Gene_Sim_Replicate[c(1:10,17:20)] %>% left_join(Simulated_Fusion_Info_2, by = c('fusion.gene.id'))
Annot_JAFFAL_3Gene_Sim_Replicate$fusionType <- mapply(function(g1, g2, g3, current_type) {
  if (current_type == "" || is.na(current_type)) {
    
    #check if reverse order tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g2, ":", g1)))) {
      return("reverse_order:tri_fusion")}
    
    #check if its a butchered tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g1, ":", g2)) |
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g3, ":", g1))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g3, ":", g2))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1, ":", g3)))) {
      return("wrong_order:tri_fusion")
    }else {
      return("false_fusion")
    }
  }
  return(current_type)}, Annot_JAFFAL_3Gene_Sim_Replicate$Gene1, Annot_JAFFAL_3Gene_Sim_Replicate$Gene2, Annot_JAFFAL_3Gene_Sim_Replicate$Gene3, Annot_JAFFAL_3Gene_Sim_Replicate$fusionType)

Annot_JAFFAL_3Gene_Sim_Replicate%>%
  mutate(depth = case_when(
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"# Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ), Replicate_Seed = case_when(
    grepl("replicate_1.", Source) ~ "1",
    grepl("replicate_2", Source) ~ "2",
    grepl("replicate_3", Source) ~ "3",
    grepl("replicate_4", Source) ~ "4",
    grepl("replicate_5", Source) ~ "5",
    grepl("replicate_6", Source) ~ "6",
    grepl("replicate_7", Source) ~ "7",
    grepl("replicate_8", Source) ~ "8",
    grepl("replicate_9", Source) ~ "9",
    grepl("replicate_10", Source) ~ "10"))

#Genion####
#Read in Data
myfiles<- list.files(path = "/bioinformatics/ryley/Algorithm_Benchmark/Replicates/Adapter_porechop_trimmed/Genion", pattern = "*_genion$", full.names = TRUE)
Genion_Sim_Replicate<- do.call(rbind, lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 0) {
    read.table(filename) %>%
      mutate(Source = basename(filename)) %>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
  } else {
    NULL
  }
}))
og_Genion <-Genion_Sim_Replicate

# Annotate
Annot_Genion_Sim_Replicate <- Genion_Sim_Replicate %>% separate(V1, into = c("V1.1", "V1.2", "V1.3"), "::", remove=FALSE) %>% separate(V2, into = c("V2.1", "V2.2", "V2.3"), "::", remove=FALSE) %>% separate(V8, into = c("chr1", "chr2", "chr3"), ";", remove=FALSE)
Annot_Genion_Sim_Replicate[c("V1", "V2")] <- lapply(Genion_Sim_Replicate[c("V1", "V2")], function(x) gsub("::", ":", x))
Annot_Genion_Sim_Replicate <- Annot_Genion_Sim_Replicate %>% left_join(Simulated_Fusion_Info_2, by = c('V1'='fusion.gene.id')) 
Annot_Genion_Sim_Replicate$fusionType <- mapply(function(g1, g2, g3, current_type, chr1, chr2, chr3, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    
    #check if this is two part or tri part fusion
    if(g3 == ""|| is.na(g3)){
      
      #check if part of tri fusion
      if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
              str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
        return("truncated_tri_fusion")} 
      else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)) |
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
        return("reverse_order:truncated_tri_fusion")} 
      
      #check for reverse order two part genes
      else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
        matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
        return(paste("reverse_order:", matching_row[1]))}
    } 
    #check for chromosomal misalignment
    else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
    }
    
    #check if this fusion contains mitochondrial genes 
    if ((grepl("chrM:", chr1, ignore.case = TRUE) & (!grepl("chrM:", chr2, ignore.case = TRUE) | !grepl("chrM:", chr3, ignore.case = TRUE))) |
        (grepl("chrM:", chr2, ignore.case = TRUE) & (!grepl("chrM:", chr1, ignore.case = TRUE)| !grepl("chrM:", chr3, ignore.case = TRUE)))|
        (grepl("chrM:", chr3, ignore.case = TRUE) & (!grepl("chrM:", chr1, ignore.case = TRUE)| !grepl("chrM:", chr2, ignore.case = TRUE)))){
      return("false_fusion:mitochondrial_genomic")} 
    
    else if (grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE) & grepl("chrM:", chr3, ignore.case = TRUE)){
      return("false_fusion:mitochondrial")}
    
    #check if its a butchered tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g1, ":", g2)) |
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g3, ":", g1))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g3, ":", g2))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1, ":", g3)))) {
      return("wrong_order:tri_fusion")}
    
    #check if reverse order tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g2, ":", g1)))) {
      return("reverse_order:tri_fusion")}
    
    #if blank and none of the above, its just a false fusion
    else if (gene_name1 == gene_name2){
      return("false_fusion:self_misalignment") 
    } else {
      return("false_fusion")}
  }
  #if it's not blank, its already been labelled
  else {return(current_type)}
}, Annot_Genion_Sim_Replicate$V1.1, Annot_Genion_Sim_Replicate$V1.2,  Annot_Genion_Sim_Replicate$V1.3, Annot_Genion_Sim_Replicate$fusionType, Annot_Genion_Sim_Replicate$chr1, Annot_Genion_Sim_Replicate$chr2 , Annot_Genion_Sim_Replicate$chr3, Annot_Genion_Sim_Replicate$V2.1, Annot_Genion_Sim_Replicate$V2.2)

Annot_Genion_Sim_Replicate <- Annot_Genion_Sim_Replicate %>%
  mutate(depth = case_when(
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))

Annot_LongGF_Sim_Replicate$Algorithm <- "LongGF"
Annot_JAFFAL_Sim$Algorithm_Replicate <- "JAFFAL"
Annot_JAFFAL_3Gene_Sim_Replicate$Algorithm <- "JAFFAL"
Annot_Genion_Sim_Replicate$Algorithm <- "Genion"
Annot_FusionSeeker_Sim_Replicate$Algorithm <- "FusionSeeker"