###############################
# Fusion Profile Stat Summary #
###############################
fusion_profile_stat_summary <-  filter(combined_data, read_supp >=2)%>%  
  dplyr::group_by(across(-read_supp)) %>%
  dplyr::summarise(read_supp = sum(read_supp), .groups = "drop") %>% 
  unique()%>% group_by(control, depth, Sequence_Identity, Algorithm, fusionType) %>%
  summarise(fusion_type_count = n(), .groups = 'drop')%>% 
  ungroup() %>% 
  complete(control, depth, Sequence_Identity, Algorithm, fusionType, fill = list(fusion_type_count = 0)) %>%
  group_by(control, depth, Sequence_Identity, Algorithm) %>%
  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall"))

fusion_profile_stat_summary$depth <- factor(fusion_profile_stat_summary$depth, levels = c("1GB", "10GB", "100GB"))
fusion_profile_stat_summary$Sequence_Identity <- factor(fusion_profile_stat_summary$Sequence_Identity, levels = c("85%", "90%", "95%"))

write_tsv(fusion_profile_stat_summary, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/fusion_profile_stat_summary.tsv")

