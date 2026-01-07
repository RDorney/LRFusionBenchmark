myfiles<-list.files(path = "/bioinformatics/siyuan/ctat_lr_results/simulated", 
                    pattern = "ctat-LR-fusion.fusion_predictions.abridged.tsv", full.names = TRUE, recursive = TRUE)
CTATLR_SIM<- bind_rows(
  lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 136) {
    read_tsv(filename,
               col_names = TRUE) %>%
      mutate(Source = basename(dirname(filename)),
             control = ifelse(grepl("Spiked", Source), "positive", "negative"))
    } else {
    NULL
  }
}))

og_CTATLR <- CTATLR_SIM

# Annotate
CTATLR_SIM_a1 <- CTATLR_SIM %>% separate(LeftBreakpoint, into = c("chrom1", "base1", "strand1"), ":", remove=FALSE) %>% separate(RightBreakpoint, into = c("chrom2", "base2", "strand2"), ":", remove=FALSE)
CTATLR_SIM_a1["#FusionName"] <- lapply(CTATLR_SIM_a1["#FusionName"], function(x) gsub("--", ":", x))
Gene_Names<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = (unique(c(CTATLR_SIM_a1$LeftGene, CTATLR_SIM_a1$RightGene))),
                       mart = ensemblv110),
                  getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                        filters = "ensembl_gene_id",
                        values = (unique(c(CTATLR_SIM_a1$LeftGene, CTATLR_SIM_a1$RightGene))),
                        mart = ensemblv110))%>% unique()
Gene_Names$chromosome_name <- paste0("chr", Gene_Names$chromosome_name)

CTATLR_SIM_a2 <- CTATLR_SIM_a1 %>% 
  left_join(Gene_Names, by= c("LeftGene" ='external_gene_name', 'chrom1'='chromosome_name')) %>%
  mutate(Gene1_ensembl_ID= coalesce(ensembl_gene_id, LeftGene)) %>% 
  select(-ensembl_gene_id)  %>% 
  left_join(Gene_Names, by= c("RightGene" ='external_gene_name', 'chrom2'='chromosome_name')) %>%
  mutate(Gene2_ensembl_ID= coalesce(ensembl_gene_id, RightGene)) %>% select(-ensembl_gene_id) %>%
  unite("fusion.gene.id" , c(Gene1_ensembl_ID, Gene2_ensembl_ID), sep = ":", remove= FALSE)
   
Annot_CTATLR_Sim <- CTATLR_SIM_a2 %>% left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id') 
#Assign label to false fusions or partially recalled fusions
Annot_CTATLR_Sim <- Annot_CTATLR_Sim
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:",first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA"))

Annot_CTATLR_Sim$fusionType <- mapply(function(g1, g2, Gen1, Gen2, current_type, chr1, chr2) {
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
}, Annot_CTATLR_Sim$Gene1_alternative_ID.x, Annot_CTATLR_Sim$Gene2_alternative_ID.x,
Annot_CTATLR_Sim$Gene1, Annot_CTATLR_Sim$Gene2, 
Annot_CTATLR_Sim$fusionType, Annot_CTATLR_Sim$Chrom1, Annot_CTATLR_Sim$Chrom2)

setdiff(unique(c(ogFusionSeeker$ID, ogFusionSeeker$Source)), unique(c(Annot_CTATLR_Sim$ID, Annot_CTATLR_Sim$Source)))

