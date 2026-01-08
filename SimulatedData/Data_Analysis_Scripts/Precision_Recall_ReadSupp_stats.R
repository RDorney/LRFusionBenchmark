#Calculate Recall vs Precision
stat_summary <- filter(combined_data, read_supp >= 2, control == "positive") %>%  unique() %>%
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
  mutate(FDR = False_Point / (False_Point + True_Point),
         Precision = True_Point / (True_Point + False_Point),
         Recall = True_Point / Spiked_Number,  
         F1 = (2 * (Precision * Recall)) / (Precision + Recall),
         Sensitivity = True_Point / (True_Point + (Spiked_Number - True_Point)))

stat_summary$depth <- factor(stat_summary$depth, levels = c("100GB", "10GB", "1GB"))
stat_summary$Sequence_Identity <- factor(stat_summary$Sequence_Identity, levels = c("85%", "90%", "95%"))

write_tsv(stat_summary, file=paste0(Analysis_output_folder, "stat_summary.tsv"))

#Calculate the effect of read support on recall/precision
calculate_read_filt_statistics <- function(data, min_read_supp) {
  data %>%
    filter(read_supp >= min_read_supp, control == "positive")  %>%  
    mutate(stat_category = case_when(
      grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType) ~ "FALSE_CALL",
      TRUE ~ "TRUE_CALL")) %>% 
    group_by(control, depth, Sequence_Identity, Algorithm, stat_category) %>%
    summarise(CALL_COUNT = n()) %>%
    summarise(False_Point = sum(CALL_COUNT[stat_category == "FALSE_CALL"], na.rm = TRUE),
              True_Point = sum(CALL_COUNT[stat_category == "TRUE_CALL"], na.rm = TRUE)) %>%
    left_join(Spiked_Counts_wo_SF[c(2:4)])%>%
    
    mutate(FDR = False_Point / (False_Point + True_Point),
           Precision = True_Point / (True_Point + False_Point),
           Recall = True_Point / Spiked_Number,  
           F1 = (2 * (Precision * Recall)) / (Precision + Recall),
           Sensitivity = True_Point / (True_Point + (Spiked_Number - True_Point)),
           minimum_read_support = min_read_supp)
}

readsupp_filtering <- calculate_read_filt_statistics(combined_data, 2) 

for(number in seq(3, max(sort(unique(combined_data$read_supp))))){
  result <- calculate_read_filt_statistics(combined_data, number)
  readsupp_filtering <- rbind(result, readsupp_filtering) 
}
readsupp_filtering$minimum_read_support <- as.numeric(readsupp_filtering$minimum_read_support)
