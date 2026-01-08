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
   
Annot_CTATLR_Sim <- CTATLR_SIM_a2 %>% 
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id', relationship = "many-to-many") 

Annot_CTATLR_Sim$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for truncated_tri_fusions
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
            str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
      return("truncated_tri_fusion")
      #check for reverse order fusions
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1))|
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
      return("reverse_order:truncated_tri_fusion") 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
      return(paste("reverse_order:", matching_row[1])) 
      #check for chromosomal misalignment
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
      #check for false fusions
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
}, Annot_CTATLR_Sim$Gene1_ensembl_ID, Annot_CTATLR_Sim$Gene2_ensembl_ID, Annot_CTATLR_Sim$fusionType, Annot_CTATLR_Sim$chrom1, Annot_CTATLR_Sim$chrom2, Annot_CTATLR_Sim$LeftGene, Annot_CTATLR_Sim$RightGene)


