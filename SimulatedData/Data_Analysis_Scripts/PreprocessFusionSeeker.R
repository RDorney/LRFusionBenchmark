#FusionSeeker_Sim # Gene1 = Gene ID 1, Gene2 = Gene ID 2
FS_DIR<-"/bioinformatics/ryley/Gencode44/Algorithms/FusionSeeker"
myfiles<-list.files(path = FS_DIR, pattern = "confident_genefusion.txt", full.names = TRUE, recursive = TRUE)
FusionSeeker_Sim <-  do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE) %>%
    mutate(Source = basename(dirname(filename)))%>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
}))
ogFusionSeeker <- do.call(rbind, lapply(myfiles, function(filename) {
  read.table(filename, header = TRUE) %>%
    mutate(Source = basename(dirname(filename)))%>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
}))

FusionSeeker_Sim <- FusionSeeker_Sim %>%
  unite("fusionGene", c(Gene1, Gene2), sep="::", remove= FALSE, na.rm = TRUE)%>%
  mutate(across(c(Gene1, Gene2), ~ gsub("\\..*","", .)))%>%
  separate(fusionGene, into= c('Gene1_vID','Gene2_vID'), sep = "::", remove = FALSE) %>%
  unite("original.fusion.gene.id", c(Gene1, Gene2), sep=":", remove= FALSE, na.rm = TRUE)
PAR_Y_problem1 <- FusionSeeker_Sim %>%
  filter(grepl("PAR_Y", Gene1_vID, ignore.case = TRUE)) %>% select(Gene1) 
PAR_Y_problem2 <-  FusionSeeker_Sim %>%
  filter(grepl("PAR_Y", Gene2_vID, ignore.case = TRUE)) %>% select(Gene2)

Gene_Name<-getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_gene_id",
                       values = (unique(c(PAR_Y_problem1$Gene1, PAR_Y_problem2$Gene2))),
                       mart = ensemblv110) %>% unique()
External_Gene_Name<-getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version",  "chromosome_name"),
                                filters = c("chromosome_name", "external_gene_name"),
                                values = list("Y", unique(Gene_Name$external_gene_name)),
                                mart = ensemblv110) %>% unique()
External_Gene_Name$chromosome_name <- paste0("chr", External_Gene_Name$chromosome_name)
FusionSeeker_Sim_ensembl <- FusionSeeker_Sim %>% left_join(select((right_join(Gene_Name, External_Gene_Name, by= "external_gene_name")), c(ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id_version, chromosome_name)), by= c('Gene1'='ensembl_gene_id.x', 'Chrom1'='chromosome_name')) %>% 
  mutate(Gene1_alternative_ID= coalesce(ensembl_gene_id.y, Gene1), Gene1_alternative_vID= coalesce(ensembl_gene_id_version, Gene1_vID)) %>%
  select(-c(ensembl_gene_id.y, ensembl_gene_id_version))%>% 
  left_join(select((right_join(Gene_Name, External_Gene_Name, by= "external_gene_name")), c(ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id_version, chromosome_name)), by= c('Gene2'='ensembl_gene_id.x', 'Chrom2'='chromosome_name')) %>% 
  mutate(Gene2_alternative_ID= coalesce(ensembl_gene_id.y, Gene2), Gene2_alternative_vID= coalesce(ensembl_gene_id_version, Gene2_vID)) %>%
  select(-c(ensembl_gene_id.y, ensembl_gene_id_version)) %>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>% unique() 

#Assign label to false fusions or partially recalled fusions
Annot_FusionSeeker_Sim <- FusionSeeker_Sim_ensembl%>%
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id')%>%
  group_by(fusionGene.x) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA"))

Annot_FusionSeeker_Sim$fusionType <- mapply(function(g1, g2, Gen1, Gen2, current_type, chr1, chr2) {
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
}, Annot_FusionSeeker_Sim$Gene1_alternative_ID.x, Annot_FusionSeeker_Sim$Gene2_alternative_ID.x,
Annot_FusionSeeker_Sim$Gene1, Annot_FusionSeeker_Sim$Gene2, 
Annot_FusionSeeker_Sim$fusionType, Annot_FusionSeeker_Sim$Chrom1, Annot_FusionSeeker_Sim$Chrom2)

setdiff(unique(c(ogFusionSeeker$ID, ogFusionSeeker$Source)), unique(c(Annot_FusionSeeker_Sim$ID, Annot_FusionSeeker_Sim$Source)))
