################################################
# Manual Filtering                             #
# Assess if using heuristics changes the stats #
################################################

calculate_read_filt_statistics_heuristic <- function(data, min_read_supp) {
  data %>%
    filter(read_supp >= min_read_supp)  %>%
    mutate(stat_category = case_when(
      grepl("false_fusion:mitochondrial|false_fusion:self_misalignment", fusionType) ~ "obvious_false",
      TRUE ~ "other")) %>% filter(stat_category == "other") %>% dplyr::select(-'stat_category') %>%  
    mutate(stat_category = case_when(
      grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType) ~ "FALSE_CALL",
      TRUE ~ "TRUE_CALL")) %>%
    
    # First summarisation: Keep fusion_type_count grouped by stat_category
    group_by(control, depth, Sequence_Identity, Algorithm, stat_category) %>%
    summarise(CALL_COUNT = n()) %>%
    summarise(False_Point = sum(CALL_COUNT[stat_category == "FALSE_CALL"], na.rm = TRUE),
              True_Point = sum(CALL_COUNT[stat_category == "TRUE_CALL"], na.rm = TRUE)) %>%
    left_join(Spiked_Counts_wo_SF[c(2:4)])%>%
    # Compute final statistics
    mutate(minimum_read_support = min_read_supp,
           FDR = False_Point / (False_Point + True_Point),
           Precision = True_Point / (True_Point + False_Point),
           Recall = True_Point / Spiked_Number,  
           F1 = (2 * (Precision * Recall)) / (Precision + Recall),
           Sensitivity = True_Point / (True_Point + (Spiked_Number - True_Point)))
}

manual_readsupp_filtering <- calculate_read_filt_statistics_heuristic(collapsed_combined_data, 2)
for(number in seq(3, max(sort(unique(collapsed_combined_data$read_supp))))){
  result <- calculate_read_filt_statistics_heuristic(collapsed_combined_data, number)
  manual_readsupp_filtering <- rbind(result, manual_readsupp_filtering) 
}
manual_readsupp_filtering$minimum_read_support <- as.numeric(manual_readsupp_filtering$minimum_read_support)
write_tsv(manual_readsupp_filtering, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/manual_readsupp_filtering.tsv")

