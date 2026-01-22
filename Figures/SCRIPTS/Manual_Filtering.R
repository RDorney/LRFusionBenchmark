################################################
# Manual Filtering                             #
# Assess if using heuristics changes the stats #
################################################
manual_filtering <- filter(combined_data, read_supp >=2) %>%  
  mutate(stat_category = case_when(
    grepl("false_fusion:mitochondrial|false_fusion:self_misalignment", fusionType) ~ "obvious_false",
    TRUE ~ "other")) %>% filter(stat_category == "other") %>% dplyr::select(-'stat_category')%>%  
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

manual_filtering$depth <- factor(manual_filtering$depth, levels = c("100GB", "10GB", "1GB"))
manual_filtering$Sequence_Identity <- factor(manual_filtering$Sequence_Identity, levels = c("85%", "90%", "95%"))

ggplot(filter(manual_filtering, control == 'positive'))+
  geom_point(aes(x = Precision, y= Recall, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  labs(title = "Recall and Precision", subtitle = "heuristic filtering", x = "Precision", y = "Recall")+
  ylim(c(0,1.00))+
  xlim(c(0,1.00))

calculate_read_filt_statistics <- function(data, min_read_supp) {
  data %>%
    filter(read_supp > min_read_supp)  %>%  
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

manual_readsupp_filtering <- calculate_read_filt_statistics(combined_data, 2) 
for(number in unique(sort(combined_data$read_supp))){
  result <- calculate_read_filt_statistics(combined_data, number)
  manual_readsupp_filtering <- rbind(result, manual_readsupp_filtering) 
}
manual_readsupp_filtering$minimum_read_support <- as.numeric(manual_readsupp_filtering$minimum_read_support)

ggplot(filter(manual_readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(title = "F1 score in relation to minimum read support", subtitle = "heuristic filtering", 
       x = "Read Support",
       y = "F1 Score")+
  ylim(c(0,1.00)) 

ggplot(filter(manual_readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(title = "Precision in relation to minimum read support", subtitle = "heuristic filtering",
       x = "Read Support",
       y = "Precision") 
