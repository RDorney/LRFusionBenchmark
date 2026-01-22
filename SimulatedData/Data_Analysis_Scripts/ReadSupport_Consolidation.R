##############################
# Read Support Consolidation #
##############################
# Fusionseeker and some of fusion callers will produce duplicate calls the same fusion multiple times with different read support because of breakpoint differences. 
# We need to consolidate these for ease in comparing how read support for the same fusion will differ by fusion caller tool
num_ALG <- length(unique(combined_data$Algorithm))
fusions_in_all_algorithms <- combined_data %>%
  filter(read_supp >=2)%>% unique() %>% #filter for read support greater than 2
  group_by(fusion.gene.id, Sequence_Identity, depth, control, fusionType) %>%
  summarise(overlap = n_distinct(Algorithm)) %>%
  filter(overlap == num_ALG) %>% #filter for fusions in all algorithms
  ungroup() %>%
  left_join(combined_data)%>% unique()%>%
  group_by(fusion.gene.id, Sequence_Identity, depth, control, Algorithm, overlap, fusionType) %>%
  summarize(Total_Read_Supp = sum(read_supp))%>%
  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall"))

fusions_in_all_algorithms$depth <- factor(fusions_in_all_algorithms$depth, levels = c("1GB", "10GB", "100GB"))
fusions_in_all_algorithms$Sequence_Identity <- factor(fusions_in_all_algorithms$Sequence_Identity, levels = c("95%", "90%", "85%"))
