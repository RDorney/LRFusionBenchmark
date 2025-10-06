#Read in fusion calling in summary data
combined_data <- read_csv("~/LongReadFusionCallerBenchmark/combined_data_FusionCalls.csv")
#Calculate Recall vs Precision
stat_summary <- filter(combined_data, read_supp > 2) %>%  
  mutate(stat_category = case_when(
    grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType) ~ "FALSE_CALL",
    TRUE ~ "TRUE_CALL")) %>%
  
  # First summarisation: Keep fusion_type_count grouped by stat_category
  group_by(control, depth, Sequence_Identity, Algorithm, stat_category) %>%
  summarise(CALL_COUNT = n()) %>%
  summarise(False_Point = sum(CALL_COUNT[stat_category == "FALSE_CALL"], na.rm = TRUE),
            True_Point = sum(CALL_COUNT[stat_category == "TRUE_CALL"], na.rm = TRUE)) %>%
  
  # Compute final statistics
  mutate(FDR = False_Point / (False_Point + True_Point),
         Precision = True_Point / (True_Point + False_Point),
         Recall = True_Point / 400,  # There are 400 spiked fusions
         F1 = (2 * (Precision * Recall)) / (Precision + Recall),
         Sensitivity = True_Point / (True_Point + (400 - True_Point)))

stat_summary$depth <- factor(stat_summary$depth, levels = c("1GB", "10GB", "100GB"))
stat_summary$Sequence_Identity <- factor(stat_summary$Sequence_Identity, levels = c("85%", "90%", "95%"))

#Calculate the effect of read support on recall/precision
calculate_read_filt_statistics <- function(data, min_read_supp) {
  data %>%
    filter(read_supp > min_read_supp)  %>%  
    mutate(stat_category = case_when(
      grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType) ~ "FALSE_CALL",
      TRUE ~ "TRUE_CALL")) %>% 
    group_by(control, depth, Sequence_Identity, Algorithm, stat_category) %>%
    summarise(CALL_COUNT = n()) %>%
    summarise(False_Point = sum(CALL_COUNT[stat_category == "FALSE_CALL"], na.rm = TRUE),
              True_Point = sum(CALL_COUNT[stat_category == "TRUE_CALL"], na.rm = TRUE)) %>%
    mutate(FDR = False_Point / (False_Point + True_Point),
           Precision = True_Point / (True_Point + False_Point),
           Recall = True_Point / 400,  # There are 400 spiked fusions
           F1 = (2 * (Precision * Recall)) / (Precision + Recall),
           Sensitivity = True_Point / (True_Point + (400 - True_Point)))
}

readsupp_filtering <- calculate_read_filt_statistics(combined_data, 2) 

for(number in unique(sort(combined_data$read_supp))){
  result <- calculate_statistics(combined_data, number)
  readsupp_filtering <- rbind(result, readsupp_filtering) 
}
readsupp_filtering$minimum_read_support <- as.numeric(readsupp_filtering$minimum_read_support)

#Generate plots
ggplot(filter(stat_summary, control == 'positive'))+
  geom_point(aes(x = Precision, y= Recall, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  labs(title = "Recall and Precision", subtitle = "minimum read support of 2", x = "Precision", y = "Recall")

ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(title = "F1 score in relation to minimum read support",
       x = "Read Support",
       y = "F1 Score") 

ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(title = "Precision in relation to minimum read support",
       x = "Read Support",
       y = "Precision") 
