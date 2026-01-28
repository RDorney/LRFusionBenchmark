JAF_DIR<-"/bioinformatics/ryley/Algorithm_Benchmark/Adapter_porechop_trimmed/JAFFAL/" 
myfiles<-list.files(path = JAF_DIR, pattern = "jaffa_results.csv", full.names = TRUE, recursive = TRUE)
JAFFAL_Sim <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = ",") %>%
    mutate(Source = basename(dirname(filename))) %>% mutate(control = ifelse(grepl("Spiked", sample), "positive", "negative"))
})) %>% separate(fusion.genes, into = c("Gene1", "Gene2"), sep = ":", remove = FALSE)
ogJAFFAL <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = ",") %>%
    mutate(Source = basename(dirname(filename))) %>% mutate(control = ifelse(grepl("Spiked", sample), "positive", "negative"))
}))

JAFFAL_Sim <- JAFFAL_Sim %>%
  filter(spanning.reads >= 2)

# fusion.genes = fusion gene name : 
Gene_Name<-getBM(attributes = c("external_gene_name", 
                                "ensembl_gene_id", 
                                "chromosome_name"),
                 filters = "external_gene_name",
                 values = unique(c(JAFFAL_Sim$Gene1, JAFFAL_Sim$Gene2)),
                 mart = ensemblv109) %>% unique()
Gene_Name$chromosome_name <- paste0("chr", Gene_Name$chromosome_name)

JAFFAL_Sim_ensembl <- left_join(JAFFAL_Sim, Gene_Name, by= c('Gene1'='external_gene_name')) %>% 
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
  dplyr::select(-c("ensembl_gene_id.x", "chromosome_name.x", "ensembl_gene_id.y", "chromosome_name.y")) %>% unique()

setdiff(unique(JAFFAL_Sim_ensembl$fusion.genes), unique(ogJAFFAL$fusion.genes))

#Annotate as known or false positive
Annot_JAFFAL_Sim <- JAFFAL_Sim_ensembl %>% 
  left_join(Simulated_Fusion_Info_2, by = c('fusion.gene.id')) %>% unique()%>%
  group_by(fusion.genes) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA")) 
# Check if ensembl_gene1 and ensembl_gene2 are contained in fusion.gene.id
Annot_JAFFAL_Sim$fusionType <- mapply(function(g1, g2, Gen1, Gen2, current_type, chr1, chr2) {
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
}, Annot_JAFFAL_Sim$ensembl_gene1, Annot_JAFFAL_Sim$ensembl_gene2, Annot_JAFFAL_Sim$Gene1, Annot_JAFFAL_Sim$Gene2, Annot_JAFFAL_Sim$fusionType, Annot_JAFFAL_Sim$chrom1, Annot_JAFFAL_Sim$chrom2)

#JAFFAL_3Gene_Sim 
myfiles<-list.files(path = JAF_DIR, pattern = ".3gene_summary", full.names = TRUE, recursive = TRUE)
JAFFAL_3Gene_Sim <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE, sep = "\t") %>%
    mutate(Source = basename(dirname(filename))) %>% 
    mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
})) %>% 
  separate(Fusion, into = c("Gene1", "Gene2", "Gene3"), sep = ":", remove = FALSE) %>%
  filter(Reads >= 2)

Gene_Name<-getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = unique(c(JAFFAL_3Gene_Sim$Gene1, JAFFAL_3Gene_Sim$Gene2, JAFFAL_3Gene_Sim$Gene3)),
                       mart = ensemblv109)

JAFFAL_3Gene_Sim <- left_join(JAFFAL_3Gene_Sim, Gene_Name, by= c('Gene1'='external_gene_name')) %>%
  left_join(Gene_Name, by= c('Gene2'='external_gene_name')) %>%
  left_join(Gene_Name, by= c('Gene3'='external_gene_name')) %>% unique() %>%
  mutate(ensembl_gene1 = ifelse(is.na(ensembl_gene_id.x) |nchar(ensembl_gene_id.x) == 0, Gene1, ensembl_gene_id.x), 
         ensembl_gene2 = ifelse(is.na(ensembl_gene_id.y) |nchar(ensembl_gene_id.y) == 0, Gene2, ensembl_gene_id.y), 
         ensembl_gene3 = ifelse(is.na(ensembl_gene_id) |nchar(ensembl_gene_id) == 0, Gene3, ensembl_gene_id)) %>%
  unite("fusion.gene.id", c(ensembl_gene1, ensembl_gene2, ensembl_gene3), sep=":", remove= FALSE, na.rm = TRUE)

Annot_JAFFAL_3Gene_Sim <- JAFFAL_3Gene_Sim[c(1:10,17:20)] %>% left_join(Simulated_Fusion_Info_2, by = c('fusion.gene.id'))
Annot_JAFFAL_3Gene_Sim$fusionType <- mapply(function(g1, g2, g3, current_type) {
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
  return(current_type)}, Annot_JAFFAL_3Gene_Sim$Gene1, Annot_JAFFAL_3Gene_Sim$Gene2, Annot_JAFFAL_3Gene_Sim$Gene3, Annot_JAFFAL_3Gene_Sim$fusionType)