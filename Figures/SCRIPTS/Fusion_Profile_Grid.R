#########################
# Fusion Profile Grid   #
#########################
# True Recall

TRC<- ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall"))+ 
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count)+
  geom_tile()+
  geom_text(aes(label = fusion_type_count), color = "black", size = 3, alpha =1) +
  facet_nested(fusionType ~ Algorithm, labeller = as_labeller(c(c("hybrid"="Inter-\n chromosomal", 
                                                                  "intra_chromosome"="Intra- \n chromosomal",
                                                                  "read_through"="Read-through",
                                                                  "tri_fusion"="Tri-fusion"), 
                                                                c("Genion"="Genion", "JAFFAL" = "JAFFAL", 
                                                                  "LongGF" = "LongGF", "FusionSeeker" = "FusionSeeker", 
                                                                  "CTAT-LR-Fusion" = "CTAT-LR-Fusion"),
                                                                c("positive"="Postive", "negative" = "Negative"))) ) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "springgreen4") +  # Adjust color scale
  scale_alpha_continuous(range = c(0.2, 1)) +  # Adjust transparency based on fusion_type_count
  labs(x = "Depth",
       y = "Mean Sequence Identity",
       fill =  "Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="lightgrey")) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/True_Recall_category.pdf", plot = TRC) #, width = 210, height = 297, units = "mm"
ggsave("~/LongReadFusionCallerBenchmark/Figures/True_Recall_category.png", plot = TRC)

# Partial Recall
PRC <- ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "Partial_Recall"))+ 
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count) +
  geom_tile() +
  geom_text(aes(label = fusion_type_count), color = "black", size = 3, alpha =1) +
  facet_nested(fusionType ~ Algorithm, labeller = as_labeller(c(c("reverse_order: hybrid"="(R) Inter-\nchromosomal", 
                                                                  "reverse_order: intra_chromosome"="(R)  Intra-\nchromosomal",
                                                                  "reverse_order: read_through"="(R) Read-through",
                                                                  "wrong_order:tri_fusion"="(W) Tri-fusion",
                                                                  "reverse_order:tri_fusion"="(R) Tri-fusion",
                                                                  "reverse_order:truncated_tri_fusion"="(R) Truncated \n Tri-fusion",
                                                                  "truncated_tri_fusion"="Truncated \n Tri-fusion" ,
                                                                  "chromosomal_misalignment:intra_chromosome"="(CM)  Intra-\nchromosomal"), 
                                                                c("Genion"="Genion", "JAFFAL" = "JAFFAL", 
                                                                  "LongGF" = "LongGF", "FusionSeeker" = "FusionSeeker", 
                                                                  "CTAT-LR-Fusion" = "CTAT-LR-Fusion"),
                                                                c("positive"="Postive", "negative" = "Negative"))) ) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "orange") +
  scale_alpha_continuous(range = c(0.2, 1)) +  # Adjust transparency based on fusion_type_count
  labs(x = "Depth",
       y = "Mean Sequence Identity",
       fill =  "#Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="lightgrey")) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/Partial_Recall_category.pdf", plot = PRC, width = 210, height = 300, units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/Partial_Recall_category.png", plot = PRC, width = 210, height = 300, units = "mm")