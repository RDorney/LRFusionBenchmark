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
Simulated_Fusion_Info_2_annotdb<-Sim_Info_Names %>%
  mutate(Gene1_alternative_ID = map(Gene1_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene2_alternative_ID = map(Gene2_vID, ~ alternative_ID[which(geneName_ID == .x)]),
         Gene3_alternative_ID = map(Gene3_vID, ~ if (!is.na(.x)) alternative_ID[which(geneName_ID == .x)] else NA_character_)) %>%
  unnest_longer(Gene1_alternative_ID) %>%
  unnest_longer(Gene2_alternative_ID) %>%
  unnest_longer(Gene3_alternative_ID, keep_empty = TRUE)%>%
  unite("fusion.gene.id" , c(Gene1_alternative_ID, Gene2_alternative_ID, Gene3_alternative_ID), sep = ":", remove= FALSE, na.rm = TRUE) %>%
  dplyr::select(-c(alternative_ID)) %>% unique() %>%
  unite("original.fusion.gene.id", c(V1, V2), sep=":", remove= FALSE, na.rm = TRUE)
Simulated_Fusion_Info_2 <- read_tsv("~/LongReadFusionCallerBenchmark/SimulatedData/simulated_data_info/simulated_fusion_info_with_GeneID.tsv")
#write_tsv(Simulated_Fusion_Info_2, file="~/LongReadFusionCallerBenchmark/SimulatedData/simulated_data_info/simulated_fusion_info_with_GeneID.tsv")

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


