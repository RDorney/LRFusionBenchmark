#R Script for producing figure 1
Figure1A <- ggplot(
  filter(combined_data_annotation, control == "positive", depth == "100GB", Sequence_Identity == "95%")
) +
  geom_bar(aes(x = Algorithm, y = True_Call_Number, fill = Algorithm), stat = "identity") +
  theme_minimal() +
  labs(title = "True Fusions",
       subtitle = "95% Sequence Identity 100Gb Depth",
       y = "Simulated Fusions",
       x = "")+
  theme(legend.position = "none")+
  ylim(c(0,400))

Figure1B <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall", depth == "100GB", Sequence_Identity == "95%"))+ 
  
  geom_bar(aes(x = Algorithm, y = fusion_type_count, fill = Algorithm), stat = "identity")+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("hybrid"="Inter-chromosomal", 
                                                     "intra_chromosome"="Intra-chromosomal",
                                                     "read_through"="Read-through",
                                                     "tri_fusion"="Tri-fusion"))))+ 
  theme_minimal() +
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")+
  ylim(c(0,100))+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Figure1AB <- (Figure1A + Figure1B)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure1AB.pdf", plot = Figure1AB)


Figure1C <- ggplot(
  filter(combined_data_annotation, control == "positive", depth == "100GB", Sequence_Identity == "95%")
) +
  geom_bar(aes(x = Algorithm, y = False_Call_Number, fill = Algorithm), stat = "identity") +
  theme_minimal() +
  labs(title = "False Fusions",
       subtitle = "95% Sequence Identity 100Gb Depth",
       y = "Fusions",
       x = "")+
  theme(legend.position = "none")

Figure1D <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "False_Call", depth == "100GB", Sequence_Identity == "95%"))+ 
  
  geom_bar(aes(x = Algorithm, y = fusion_type_count, fill = Algorithm), stat = "identity")+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order Tri-fusion", 
                                                     "false_fusion:self_misalignment" = "Self-Misalignment",
                                                     "false_fusion:mitochondrial_genomic" = "Mitochondrial:Genomic",
                                                     "false_fusion:mitochondrial" = "Mitochondrial",
                                                     "false_fusion" = "False Chimera"))))+ 
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")

Figure1CD <- (Figure1C + Figure1D)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure1CD.pdf", plot = Figure1CD)

#Overlap true fusions plot####
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
Figure1E <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:4]), name='',
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='white'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='springgreen4'), guide='none')
    + ylab('True Recall')
    ))) 
ggsave("Compare_Contrast_True_Calls.pdf", plot = Figure1E, width = 297, height = 210, units = "mm")

Figure1 <- Figure1AB | Figure1CD)/
  (Figure1E) +   
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure1.pdf", plot = Figure1)