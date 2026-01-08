###############################
# Fusion Profile Stat Summary #
###############################
fusion_profile_stat_summary <- combined_data %>%  filter(read_supp >=2) %>% unique()%>% group_by(control, depth, Sequence_Identity, Algorithm, fusionType) %>%
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

FCC <- ggplot(filter(fusion_profile_stat_summary, recall_category == "False_Call"))+ #, control == "positive"
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count) +
  geom_tile() +
  geom_text(aes(label = fusion_type_count), color = "black", size = 3, alpha =1) +
  facet_nested(control + fusionType ~ Algorithm, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order Tri-fusion", 
                                                                            "false_fusion:self_misalignment" = "Self-Misalignment",
                                                                            "false_fusion:mitochondrial_genomic" = "Mitochondrial:Genomic",
                                                                            "false_fusion:mitochondrial" = "Mitochondrial",
                                                                            "false_fusion" = "False Chimera"), 
                                                                          c("Genion"="Genion", "JAFFAL" = "JAFFAL", "LongGF" = "LongGF", "FusionSeeker" = "FusionSeeker"),
                                                                          c("positive"="Postive", "negative" = "Negative"))) ) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "red")+
  scale_alpha_continuous(range = c(0.2, 1)) +  # Adjust transparency based on fusion_type_count
  labs(x = "Depth",
       y = "Mean Sequence Identity",
       fill =  "#Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="lightgrey")) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/False_fusion_category_counts.pdf", plot = FCC, width = 210, height = 297, units = "mm")