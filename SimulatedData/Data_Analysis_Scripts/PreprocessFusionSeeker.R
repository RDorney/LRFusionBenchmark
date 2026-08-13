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
  unite("fusion.gene.id", c(Gene1, Gene2), sep=":", remove= FALSE, na.rm = TRUE) %>%
  filter(NumSupp >= 2)
#########################################
# Check FusionSeeker for antisense genes
#########################################
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
# get unique genes from both columns
genes_to_check <- unique(c(FusionSeeker_Sim$Gene1, FusionSeeker_Sim$Gene2))

# run get_antisense_overlaps on each gene
FusionSeeker_antisense_list <- lapply(genes_to_check, 
                                      function(g) get_antisense_overlaps(g, genes, id_type = "gene_id"))

# combine into a single data frame
FusionSeeker_antisense_df <- bind_rows(FusionSeeker_antisense_list, .id = "query_gene_index")

# add the original gene names
FusionSeeker_antisense_df$query_gene <- genes_to_check[as.integer(FusionSeeker_antisense_df$query_gene_index)]
FusionSeeker_antisense_df$query_gene_index <- NULL

antisense_pairs <- FusionSeeker_antisense_df %>%
  dplyr::mutate(target_gene_id = sub("\\..*$", "", target_gene_id),
                gene_id        = sub("\\..*$", "", gene_id)) %>%
  dplyr::select(target_gene_id, gene_id) %>%
  distinct()

#############################################################
#Assign label to false fusions or partially recalled fusions
#############################################################
Annot_FusionSeeker_Sim <- FusionSeeker_Sim %>%
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id')%>%
  group_by(fusionGene.x) %>%
  mutate(fusionType = case_when(
    is.na(fusionType) ~ paste0("chromosomal_misalignment:", dplyr::first(na.omit(fusionType))),  # Only update NA values
    .default = fusionType  # Keep existing values unchanged
  ))%>% ungroup()%>%
  replace_with_na(replace = list(fusionType = "chromosomal_misalignment:NA"))

Annot_FusionSeeker_Sim$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
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
      return(paste0("reverse_order:", matching_row[1])) 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste0("chromosomal_misalignment:", matching_row[1])) 
    } else if ((grepl("chrM", chr1, ignore.case = TRUE) & !grepl("chrM", chr2, ignore.case = TRUE)) | (!grepl("chrM", chr1, ignore.case = TRUE) & grepl("chrM", chr2, ignore.case = TRUE))){
      return("false_fusion:mitochondrial_genomic") 
    }else if (grepl("chrM", chr1, ignore.case = TRUE) & grepl("chrM", chr2, ignore.case = TRUE)){
      return("false_fusion:mitochondrial") 
    }else if (((g1 == g2)) & (chr1 != chr2)){
      return("false_fusion:self_misalignment") 
    }else if (paste(g1, g2) %in% paste(antisense_pairs$target_gene_id, antisense_pairs$gene_id)
              | 
              paste(g2, g1) %in% paste(antisense_pairs$target_gene_id, antisense_pairs$gene_id)){
      return("false_fusion:Sense-Antisense") 
    }else {
      return("false_fusion")
    }
  }
  return(current_type)
}, Annot_FusionSeeker_Sim$Gene1, Annot_FusionSeeker_Sim$Gene2,
Annot_FusionSeeker_Sim$fusionType, Annot_FusionSeeker_Sim$Chrom1, Annot_FusionSeeker_Sim$Chrom2)

setdiff(unique(c(ogFusionSeeker$ID, ogFusionSeeker$Source)), unique(c(Annot_FusionSeeker_Sim$ID, Annot_FusionSeeker_Sim$Source)))
