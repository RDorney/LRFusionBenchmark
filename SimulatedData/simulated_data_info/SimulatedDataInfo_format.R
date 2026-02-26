simulated_fusions <- read.table(file = "/home/ryleyd/LongReadFusionCallerBenchmark/SimulatedData/simulated_data_info/FUSIM_benchmark_fusions.txt", header=TRUE, sep= "\t")
#convert from txt to tsv and save file
write_tsv(simulated_fusions, file="~/simulated_fusion_info.tsv")
# Define the standard set (1-22, X, Y, MT)
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "MT", "M"))
#seperate the fusion gene IDs into individual gene IDs for downstream annotation steps
Simulated_Fusion_Info <- simulated_fusions %>% 
  separate(fusionGene, into= c('V1','V2','V3'), sep = "-", remove = FALSE) %>% 
  mutate(across(fusionGene, ~ gsub("-", "::", .))) %>% 
  mutate(geneName_ID = geneName) %>% 
  mutate(across(c(V1, V2, V3, geneName), ~ gsub("\\..*","", .)))%>% 
  separate(fusionGene, into= c('Gene1_vID','Gene2_vID','Gene3_vID'), sep = "::", remove = FALSE)
# Retrieve gene names for the corresponding Ensembl IDs. 
# These mappings are essential for accurate comparison between simulated and experimentally detected fusions across different annotation systems.
ah <- AnnotationHub()
edb109 <- query(ah, c("EnsDb", "Homo sapiens", "109"))[[1]]
edb110 <- query(ah, c("EnsDb", "Homo sapiens", "110"))[[1]]

map_genes <- function(edb) {
  genes(edb, return.type = "data.frame") %>%
    dplyr::select(gene_name, gene_id, gene_id_version, seq_name) %>%
    mutate(seq_name = paste0("chr", seq_name))
}
combined_ref <- rbind(map_genes(edb109), map_genes(edb110))%>%
  filter(seq_name %in% standard_chrs) %>% unique()

Sim_Info_Names <- Simulated_Fusion_Info %>%
  # Match original geneName to Ensembl ID to get external name
  left_join(combined_ref %>% dplyr::select(gene_name, gene_id), 
            by = c("geneName" = "gene_id")) %>%
  # Match by name and chromosome to find alternative IDs
  left_join(combined_ref, 
            by = c("gene_name", "chrom" = "seq_name"), 
            suffix = c("", ".alt")) %>%
  mutate(alternative_ID = coalesce(gene_id, geneName)) %>%
  dplyr::select(-gene_id) %>%
  rename(external_gene_name = gene_name)

#create a dataframe inclusive of alternative fusion gene names as gene names can be variable.  
Simulated_Fusion_Info_2<-Sim_Info_Names %>%
  mutate(Gene1_alternative_ID = map(Gene1_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene2_alternative_ID = map(Gene2_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene3_alternative_ID = map(Gene3_vID, ~ if (!is.na(.x)) alternative_ID[which(geneName_ID == .x)] else NA_character_)) %>%
  unnest_longer(Gene1_alternative_ID) %>%
  unnest_longer(Gene2_alternative_ID) %>%
  unnest_longer(Gene3_alternative_ID, keep_empty = TRUE)%>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID, Gene3_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>%
  dplyr::select(-c(alternative_ID)) %>% unique() %>%
  unite("original.fusion.gene.id", c(V1, V2), sep=":", remove= FALSE, na.rm = TRUE)
write_tsv(Simulated_Fusion_Info_2, file="~/LongReadFusionCallerBenchmark/SimulatedData/simulated_data_info/simulated_fusion_info_with_GeneID.tsv")

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

#convert from txt to tsv and save file
write_tsv(Simulated_Fusion_Info_2, file="~/LongReadFusionCallerBenchmark/SimulatedData/simulated_data_info/simulated_fusion_info_with_GeneID.tsv")


#################################
# Old method with BioMaRt
#################################
Gene1_Name<-rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
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
                                mart = ensemblv109), 
                          getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name"),
                                filters = "external_gene_name",
                                values = unique(Gene1_Name$external_gene_name),
                                mart = ensemblv110)) %>% unique()
External_Gene_Name$chromosome_name <- paste0("chr",External_Gene_Name$chromosome_name)

Alternative_Names_Gene_Fusion <- rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
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
    alternative_ID = coalesce(ensembl_gene_id, geneName)) %>% dplyr::select(-c(ensembl_gene_id))

