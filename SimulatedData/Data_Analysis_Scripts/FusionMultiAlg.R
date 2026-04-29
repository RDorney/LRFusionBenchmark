##################################
# Fusions in Multiple Algorithms #
##################################
# identify fusions in more than one algorithm
#consolidate reads for duplicate fusion calls due to different breakpoints
Total_Read_Supp_summary <- filter(combined_data, read_supp >= 2) %>%  
  dplyr::group_by(across(-read_supp)) %>%
  dplyr::summarise(read_supp = sum(read_supp), .groups = "drop")%>% unique() %>% #remove duplicates as a result of annotations
  group_by(fusion.gene.id, Sequence_Identity, depth, control, Algorithm, fusionType)  %>%
  summarise(Total_Read_Supp = sum(read_supp), .groups = "drop")

# Second summarise
Algorithm_summary <- Total_Read_Supp_summary %>%
  group_by(fusion.gene.id, Sequence_Identity, depth, control, fusionType) %>%
  summarise(overlap = n_distinct(Algorithm), .groups = "drop")

# Join to retain Total_Read_Supp
fusions_in_n_algorithms <- left_join(Total_Read_Supp_summary, Algorithm_summary, 
                                     by = c("fusion.gene.id", "Sequence_Identity", "depth", "control", "fusionType")) %>%
  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall"))
fusions_in_n_algorithms$depth <- factor(fusions_in_n_algorithms$depth, levels = c("1GB", "10GB", "100GB"))
fusions_in_n_algorithms$Sequence_Identity <- factor(fusions_in_n_algorithms$Sequence_Identity, levels = c("95%", "90%", "85%"))
fusions_in_n_algorithms <- fusions_in_n_algorithms %>%
  mutate(fusion.gene.id = fct_reorder(fusion.gene.id, overlap, .desc = TRUE))
write_tsv(fusions_in_n_algorithms, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/fusions_in_n_algorithms.tsv")
