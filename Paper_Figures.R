##### Library
library(patchwork)
library(ggplot2)
library(ggforce)
library(ComplexUpset)
#Figure 2 ####
#Fusion Abundance Plots
Figure2A_1 <- ggplot(
  filter(combined_data_annotation, control == "positive", depth == "100GB", Sequence_Identity == "95%")
) +
  geom_bar(aes(x = Algorithm, y = True_Call_Number, fill = Algorithm), stat = "identity") +
  theme_minimal() +
  labs(title = "True Fusions",
       subtitle = "95% Sequence Identity 100Gb Depth",
       y = "Simulated Fusions",
       x = "")+
  theme(legend.position = "none")+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))+
  ylim(c(0,400))

Figure2A_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall", depth == "100GB", Sequence_Identity == "95%"))+ 
  
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
        axis.ticks.x = element_blank())+
theme(
  plot.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 12))+

Figure2A <- (Figure2A_1 + Figure2A_2)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure2A.pdf", plot = Figure2A)


Figure2B_1 <- ggplot(
  filter(combined_data_annotation, control == "positive", depth == "100GB", Sequence_Identity == "95%")
) +
  geom_bar(aes(x = Algorithm, y = False_Call_Number, fill = Algorithm), stat = "identity") +
  theme_minimal() +
  labs(title = "False Chimeras",
       subtitle = "95% Sequence Identity 100Gb Depth",
       y = "Chimeras",
       x = "")+
  theme(legend.position = "none")+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

Figure2B_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "False_Call", depth == "100GB", Sequence_Identity == "95%"))+ 
  
  geom_bar(aes(x = Algorithm, y = fusion_type_count, fill = Algorithm), stat = "identity")+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order Tri-fusion", 
                                                     "false_fusion:self_misalignment" = "Self\nMisalignment",
                                                     "false_fusion:mitochondrial_genomic" = "Mitochondrial\nGenomic",
                                                     "false_fusion:mitochondrial" = "Mitochondrial",
                                                     "false_fusion" = "False Chimera"))))+ 
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")

Figure2B <- (Figure2B_1 + Figure2B_2)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure2B.pdf", plot = Figure2B)

#Overlap true fusions plots####
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp > 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
True_overlap <- ComplexUpset::upset(
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
ggsave("Compare_Contrast_True_Calls.pdf", plot = True_overlap, width = 297, height = 210, units = "mm")

#Figure2s####
Figure2_paper <- (Figure2A | Figure2B)/
                  (True_overlap) +   
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure2_paper.pdf", plot = Figure2_paper)

#Figure 3 ####
stat_summary_2<- stat_summary
stat_summary_2$depth <- factor(stat_summary_2$depth, levels = c("1GB", "10GB", "100GB"))

Figure3A_1 <- ggplot(filter(
  stat_summary_2, control == "positive", Sequence_Identity == "95%"),
  aes(x = depth, y = Precision, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity))) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequencing Depth",
       subtitle = "95% Sequence Identity",
       y = "Precision",
       x = "")+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
Figure3A_2 <-ggplot(filter(
  stat_summary_2, control == "positive", depth == "100GB"),
  aes(x = Sequence_Identity, y = Precision, colour = Algorithm, group = interaction(Algorithm, depth))) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequence Identity",
       subtitle = "100Gb Sequencing Depth",
       y = "",
       x = "")+
  ylim(c(0,1.00))+
  theme(legend.position = "none",
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
Figure3B_1 <-ggplot(filter(
  stat_summary_2, control == "positive", Sequence_Identity == "95%"),
  aes(x = depth, y = Recall, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity))) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "",
       subtitle = "",
       y = "Recall",
       x = "Depth")+
  ylim(c(0,1.00))+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

Figure3B_2 <-ggplot(filter(
  stat_summary_2, control == "positive", depth == "100GB"),
  aes(x = Sequence_Identity, y = Recall, colour = Algorithm, group = interaction(Algorithm, depth))) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "",
       subtitle = "",
       y = "",
       x = "Sequence Identity")+
  ylim(c(0,1.00))+
  theme(legend.position = "none",
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

Figure3 <- (Figure3A_1 + Figure3A_2)/
  (Figure3B_1 + Figure3B_2)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("Figure3_paper.pdf", plot = Figure3, width = 297, height = 210, units = "mm")

#Figure 4 ####
F1_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "F1")+
  ylim(c(0,1.00)) +
  theme(legend.position = "right", 
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
Precision_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "", y = "Precision") +
  ylim(c(0,1.00))+
  theme(legend.position = "none",
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)) +
  labs(title = "Sequencing Depth 100Gb , 95% Mean Sequence Identity")
Recall_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = Recall, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "Recall")+
  ylim(c(0,1.00))+
  theme(legend.position = "none",
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

Figure4A <-((Precision_subset /  Recall_subset + plot_layout(guides = 'keep')) |  F1_subset)
Figure4B <-ggplot(filter(fusions_in_n_algorithms, depth == "100GB", Sequence_Identity == "95%", control == "positive", overlap == 4), aes(x = fusion.gene.id, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 5, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="",
       fill = "Read Support") +
  theme_minimal() +
  facet_nested(Sequence_Identity ~ depth + recall_category , scales = "free_x", space = "free")+
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 50)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = "left",
        legend.spacing.x = unit(0.1, "cm"),   # reduce horizontal gap
        legend.spacing.y = unit(0.1, "cm"),   # reduce vertical gap inside legend
        legend.margin    = margin(0, 0, 0, 0),# trim around the legend box
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))+ plot_layout(guides = 'keep')
Figure4A + Figure4B
ggsave("Figure4A_paper.pdf", plot = Figure4A, width = 297, height = 210, units = "mm")
ggsave("Figure4B_paper.pdf", plot = Figure4B, width = 297, height = 210, units = "mm")  

#Figure 5 ####
Figure5 <- ggplot(filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%"), aes(y = values_vector, x = Breakpoint_number, color = Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina(maxwidth = 0.8)+
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),   # remove minor log lines
    panel.grid.major.x = element_blank(),  # optional: remove vertical major lines
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))+ 
  scale_y_continuous(trans='log1p', labels = label_log(digits = 2),
                     breaks = c(0, 1, 10, 100, 1000, 10000, 100000, 1000000, 1000000, 10000000, 100000000))+
  labs(y = "Absolute distance from simulated breakpoint",
       x= "")

ggsave("Figure5_paper.pdf", plot = Figure5, width = 297, height = 210, units = "mm")

#Figure 6 ####
# Plot for K562 with smaller y-axis limit
p_k562 <- ggplot(filter(super_summary, Sequencing_Depth != "Total", Cell_Lines == "K562")) +
  geom_point(aes(x = novelnumber, y = knownnumber, colour = Algorithm, shape = Sequencing_Depth)) +
  facet_grid(Library ~ Cell_Lines)+ 
  theme(legend.position = "none",
        strip.text.y = element_blank(),         # removes row facet labels (Library)
        strip.background.y = element_blank(),    # removes grey background for row labels
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(x = "Novel Fusions", y = "Known Fusions") +
  ylim(0, 6)  # Adjust this limit to your preference

# Plot for MCF7 with larger y-axis limit
p_mcf7 <- ggplot(filter(super_summary, Sequencing_Depth != "Total", Cell_Lines == "MCF7")) +
  geom_point(aes(x = novelnumber, y = knownnumber, colour = Algorithm, shape = Sequencing_Depth)) +
  facet_grid(Library ~ Cell_Lines) +
  labs(x = "Novel Fusions", y = "")+ 
  theme(strip.text.y = element_blank(),         # removes row facet labels (Library)
        strip.background.y = element_blank(),    # removes grey background for row labels
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  ylim(0, 53)  # Adjust as needed

#combine the two plots
Figure6 <- p_k562 + p_mcf7
ggsave("Figure6_paper.pdf", plot = Figure6, width = 210, height =200 , units = "mm")
ggsave("Figure6_paper.jpeg", plot = Figure6, width = 210, height =200 , units = "mm")

#Supplementary Figures 1 & 2  ####
#supp fig 1
TRC<- ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall"))+ 
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count)+
  geom_tile()+
  geom_text(aes(label = fusion_type_count), color = "black", size = 3, alpha =1) +
  facet_nested(fusionType ~ Algorithm, labeller = as_labeller(c(c("hybrid"="Inter-chromosomal", 
                                                                  "intra_chromosome"="Intra-chromosomal",
                                                                  "read_through"="Read-through",
                                                                  "tri_fusion"="Tri-fusion"), 
                                                                c("Genion"="Genion", "JAFFAL" = "JAFFAL", "LongGF" = "LongGF", "FusionSeeker" = "FusionSeeker"),
                                                                c("positive"="Postive", "negative" = "Negative"))) ) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "springgreen4") +  # Adjust color scale
  scale_alpha_continuous(range = c(0.2, 1)) +  # Adjust transparency based on fusion_type_count
  labs(x = "Depth",
       y = "Mean Sequence Identity",
       fill =  "#Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="lightgrey")) 
ggsave("True_Recall_category.pdf", plot = TRC, width = 210, height = 297, units = "mm")

#sup fig 2
combined_data_annotation <- filter(combined_data, read_supp > 2) %>%  mutate(recall_category = case_when(
  grepl("false", fusionType) ~ "False_Call",
  grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
  TRUE ~ "True_Recall"))%>%
  group_by(control, depth, Sequence_Identity, Algorithm, recall_category) %>%
  summarise(RECALL_COUNT = n()) %>%
  summarise(False_Call_Number = sum(RECALL_COUNT[recall_category == "False_Call"], na.rm = TRUE),
            Partial_Call_Number = sum(RECALL_COUNT[recall_category == "Partial_Recall"], na.rm = TRUE),
            True_Call_Number = sum(RECALL_COUNT[recall_category == "True_Recall"], na.rm = TRUE))%>% 
  ungroup() %>% 
  complete(control, depth, Sequence_Identity, Algorithm, fill = list(False_Call_Number = 0, Partial_Call_Number = 0, True_Call_Number = 0))
combined_data_annotation$depth <- factor(combined_data_annotation$depth, levels = c("1GB", "10GB", "100GB"))
combined_data_annotation$Sequence_Identity <- factor(combined_data_annotation$Sequence_Identity, levels = c("85%", "90%", "95%"))

write_tsv(combined_data_annotation, file="/bioinformatics/ryley/Algorithm_Benchmark/benchmark_data_analysis/Analysis/benchmark_fusion_results/combined_data_annotation.tsv")

false_call_plot <-ggplot(filter(combined_data_annotation), aes(x = depth, y = Sequence_Identity, fill = False_Call_Number)) +
  geom_tile() +
  geom_text(aes(label = False_Call_Number), color = "black", size = 3) +
  facet_grid( control ~ Algorithm) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "red")+  # Adjust color scale
  labs(x = "",
       y = "",
       fill = "False Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="grey85"))

partial_call_plot <- ggplot(filter(combined_data_annotation, control =="positive"), aes(x = depth, y = Sequence_Identity, fill = Partial_Call_Number)) +
  geom_tile() +
  geom_text(aes(label = Partial_Call_Number), color = "black", size = 3) +
  facet_grid(control ~ Algorithm)+
  scale_fill_gradient(low = "white", high = "orange") +  # Adjust color scale
  labs(x = "",
       y = "Mean Sequence Identity",
       fill = "Partially Recalled Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="grey85"))

true_call_plot <-ggplot(filter(combined_data_annotation, control =="positive"), aes(x = depth, y = Sequence_Identity, fill = True_Call_Number)) +
  geom_tile() +
  geom_text(aes(label = True_Call_Number), color = "black", size = 3) +
  facet_grid(control ~ Algorithm) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "springgreen4") +  # Adjust color scale
  labs(x = "Sequencing Depth",
       y = "",
       fill = "True Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="grey85"))

final_plot <- (false_call_plot / partial_call_plot / true_call_plot)+ 
  plot_layout(heights = c(2, 1, 1)) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))

# Display the combined plot
final_plot

#removed sup fig3
FCC <- ggplot(filter(fusion_profile_stat_summary, recall_category == "False_Call"))+ #, control == "positive"
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count) +
  geom_tile() +
  geom_text(aes(label = fusion_type_count), color = "black", size = 3, alpha =1) +
  facet_nested(control + fusionType ~ Algorithm, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order \n Tri-fusion", 
                                                                            "false_fusion:self_misalignment" = "Self\nMisalignment",
                                                                            "false_fusion:mitochondrial_genomic" = "Mitochondrial\nGenomic",
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
ggsave("False_fusion_category_counts.pdf", plot = FCC, width = 210, height = 297, units = "mm")
#Supplementary figure 3####
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp > 2)$fusion.gene.id),
  Genion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp > 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp > 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp > 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
# Generate an upset plot without specifying intersect, which uses the default behavior
PartialCall_Overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:4]),
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='orange'), guide='none')
    + ylab('Partial Recall')
    ))) 
#False overlap
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("false", fusionType), read_supp > 2)$fusion.gene.id),
  Genion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("false", fusionType), read_supp > 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("false", fusionType), read_supp > 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("false", fusionType), read_supp > 2)$fusion.gene.id)
)

Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)

# Convert list to binary membership matrix
Fusion_overlap_matrix <- as.data.frame(Fusion_overlap)  # Convert list-based to a dataframe

# Ensure the matrix is in the correct format (0/1 values)
Fusion_overlap_matrix <- as.data.frame(lapply(Fusion_overlap_matrix, as.integer))

# Generate the UpSet plot using ComplexUpset
False_Overlap <- ComplexUpset::upset(
  Fusion_overlap_matrix, names(Fusion_overlap_matrix)[1:4],  # First four sets
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='red'), guide='none') 
    + ylab('False Calls')
    )))

supp_Figure3 <-(False_Overlap/PartialCall_Overlap)
ggsave("supp_Figure3.pdf", plot = supp_Figure3, width = 297, height =210 , units = "mm")

#Supplementary figure 4####
PRC <- ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "Partial_Recall"))+ 
  aes(x = depth, y = Sequence_Identity, fill = fusion_type_count) +
  geom_tile() +
  geom_text(aes(label = fusion_type_count), color = "black", size = 4, alpha =1) +
  facet_nested(fusionType ~ Algorithm, labeller = as_labeller(c(c("reverse_order: hybrid"="(R) Inter-\nchromosomal", 
                                                                  "reverse_order: intra_chromosome"="(R) Intra-\nchromosomal",
                                                                  "reverse_order: read_through"="(R) Read-through",
                                                                  "wrong_order:tri_fusion"="(W) Tri-fusion",
                                                                  "reverse_order:tri_fusion"="(R) Tri-fusion",
                                                                  "reverse_order:truncated_tri_fusion"="(R) Truncated \nTri-fusion",
                                                                  "truncated_tri_fusion"="Truncated \n Tri-fusion" ,
                                                                  "chromosomal_misalignment:intra_chromosome"="(CM) Intra-\nchromosomal"), 
                                                                c("Genion"="Genion", "JAFFAL" = "JAFFAL", "LongGF" = "LongGF", "FusionSeeker" = "FusionSeeker"),
                                                                c("positive"="Postive", "negative" = "Negative"))) ) +  # Separate heatmaps for each algorithm
  scale_fill_gradient(low = "white", high = "orange") +
  scale_alpha_continuous(range = c(0.2, 1)) +  # Adjust transparency based on fusion_type_count
  labs(x = "Depth",
       y = "Mean Sequence Identity",
       fill =  "Fusions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.background = element_rect(color="black", fill="lightgrey"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) 
ggsave("Partial_Recall_category.pdf", plot = PRC, width = 210, height = 297, units = "mm")
#Supplementary figure 5####
supp_Figure5A <- ggplot(
  filter(combined_data_annotation, control == "positive", depth == "100GB", Sequence_Identity == "95%")
) +
  geom_bar(aes(x = Algorithm, y = Partial_Call_Number, fill = Algorithm), stat = "identity") +
  theme_minimal() +
  labs(title = "Partially Recalled Fusions",
       subtitle = "95% Sequence Identity 100Gb Depth",
       y = "Simulated Fusions",
       x = "")+
  theme(legend.position = "none")+
  ylim(c(0,400))+
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

supp_Figure5B <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "Partial_Recall", depth == "100GB", Sequence_Identity == "95%"))+ 
  
  geom_bar(aes(x = Algorithm, y = fusion_type_count, fill = Algorithm), stat = "identity")+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("reverse_order: hybrid"="(R) Inter-Chromosomal", 
                                                     "reverse_order: intra_chromosome"="(R) Intra-Chromosomal",
                                                     "reverse_order: read_through"="(R) Read-through",
                                                     "wrong_order:tri_fusion"="(W) Tri-Fusion",
                                                     "reverse_order:tri_fusion"="(R) Tri-Fusion",
                                                     "reverse_order:truncated_tri_fusion"="(R) Truncated Tri-Fusion",
                                                     "truncated_tri_fusion"="Truncated Tri-Fusion" ,
                                                     "chromosomal_misalignment:intra_chromosome"="(CM) Intra-Chromosomal"))))+ 
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_blank())+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")+
  ylim(c(0,100))+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

supp_Figure5 <- (supp_Figure5A + supp_Figure5B)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")), tag_levels = 'A')
ggsave("supp_Figure5.pdf", plot = supp_Figure5, width = 297, height = 210, units = "mm")
ggsave("supp_Figure5.png", plot = supp_Figure5, width = 297, height = 210, units = "mm")

#Supplementary figure 6####
supp_Figure6A_1 <-ggplot(filter(
  combined_data_annotation, control == "positive", Sequence_Identity == "95%"),
  aes(x = depth, y = True_Call_Number, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity))
) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequencing Depth ~ True Fusions",
       subtitle = "95% Sequence Identity",
       y = "Simulated Fusions",
       x = "")+
  ylim(c(0,400))+ 
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

supp_Figure6A_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall", Sequence_Identity == "95%"),
                    aes(x = depth, y = fusion_type_count, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity)))+ 
  geom_point() +
  geom_line()+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("hybrid"="Inter-chromosomal", 
                                                     "intra_chromosome"="Intra-chromosomal",
                                                     "read_through"="Read-through",
                                                     "tri_fusion"="Tri-fusion"))))+ 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")+
  ylim(c(0,100))

supp_Figure6B_1 <-ggplot(filter(
  combined_data_annotation, control == "positive", Sequence_Identity == "95%"),
  aes(x = depth, y = False_Call_Number, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity))
) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequencing Depth ~ False Fusions",
       subtitle = "95% Sequence Identity",
       y = "Fusions",
       x = "")+ 
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

supp_Figure6B_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "False_Call", Sequence_Identity == "95%"),
                    aes(x = depth, y = fusion_type_count, colour = Algorithm, group = interaction(Algorithm, Sequence_Identity)))+ 
  geom_point() +
  geom_line()+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order Tri-fusion", 
                                                     "false_fusion:self_misalignment" = "Self-Misalignment",
                                                     "false_fusion:mitochondrial_genomic" = "Mitochondrial:Genomic",
                                                     "false_fusion:mitochondrial" = "Mitochondrial",
                                                     "false_fusion" = "False Chimera"))))+ 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")

supp_Figure6 <- ((supp_Figure6A_1 + supp_Figure6A_2))/
  ((supp_Figure6B_1 + supp_Figure6B_2))+ plot_layout(tag_level = 'new')+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("supp_Figure6.pdf", plot = supp_Figure6,  width = 297, height = 210 , unit ="mm")

#Supplementary figure 7####
supp_Figure7A_1 <-ggplot(filter(
  combined_data_annotation, control == "positive", depth == "100GB"),
  aes(x = Sequence_Identity, y = True_Call_Number, colour = Algorithm, group = interaction(Algorithm, depth))
) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequence Identity ~ True Fusions",
       subtitle = "100Gb Sequencing Depth",
       y = "Simulated Fusions",
       x = "")+
  ylim(c(0,400))

supp_Figure7A_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "True_Recall", depth == "100GB"),
                    aes(x = Sequence_Identity, y = fusion_type_count, colour = Algorithm, group = interaction(Algorithm, depth)))+ 
  geom_point() +
  geom_line()+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("hybrid"="Inter-chromosomal", 
                                                     "intra_chromosome"="Intra-chromosomal",
                                                     "read_through"="Read-through",
                                                     "tri_fusion"="Tri-fusion"))))+ 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")+
  ylim(c(0,100))+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

supp_Figure7B_1 <-ggplot(filter(
  combined_data_annotation, control == "positive", depth == "100GB"),
  aes(x = Sequence_Identity, y = False_Call_Number, colour = Algorithm, group = interaction(Algorithm, depth))
) +
  geom_point() +
  geom_line()+
  theme_minimal() +
  labs(title = "Sequence Identity ~ False Fusions",
       subtitle = "100Gb Sequencing Depth",
       y = "Fusions",
       x = "")+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

supp_Figure7B_2 <-ggplot(filter(fusion_profile_stat_summary, control == "positive", recall_category == "False_Call", depth == "100GB"),
                    aes(x = Sequence_Identity, y = fusion_type_count, colour = Algorithm, group = interaction(Algorithm, depth)))+ 
  geom_point() +
  geom_line()+
  facet_wrap(~fusionType, labeller = as_labeller(c(c("wrong_order:tri_fusion" = "Wrong Order Tri-fusion", 
                                                     "false_fusion:self_misalignment" = "Self-Misalignment",
                                                     "false_fusion:mitochondrial_genomic" = "Mitochondrial:Genomic",
                                                     "false_fusion:mitochondrial" = "Mitochondrial",
                                                     "false_fusion" = "False Chimera"))))+ 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(title = "",
       subtitle = "",
       y = "",
       x = "")

supp_Figure7 <- (supp_Figure7A_1 +supp_Figure7A_2)/
  (supp_Figure7B_1 + supp_Figure7B_2)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave("supp_Figure7.pdf", plot = supp_Figure7, width = 297, height = 210 , unit ="mm")

#Supplementary figure 8####
Recall_Precision <- ggplot(filter(stat_summary, control == 'positive'))+
  geom_point(aes(x = Precision, y= Recall, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  labs(x = "Precision", y = "Recall")+
  ylim(c(0,1.00))+
  xlim(c(0,1.00)) +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("Recall-Precision.pdf", plot = Recall_Precision, width = 297, height = 210, units = "mm")
ggsave("Recall-Precision.png", plot = Recall_Precision, width = 297, height = 210, units = "mm")
#Supplementary figure 9####
ReadSupport_Precision <- ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "Precision") +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("ReadSupport-Precision.pdf", plot = ReadSupport_Precision, width = 297, height = 210, units = "mm")
ggsave("ReadSupport-Precision.png", plot = ReadSupport_Precision, width = 297, height = 210, units = "mm")
#Supplementary figure 10 ####
ReadSupport_F1Score <- ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "F1 Score")+
  ylim(c(0,1.00))+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) 
ggsave("ReadSupport-F1Score.pdf", plot = ReadSupport_F1Score, width = 297, height = 210, units = "mm")
#Supplementary figure 11 ####
supp_Figure11_A<- ggplot(filter(fusions_in_n_algorithms, depth == "100GB", control == "positive", overlap == 4), aes(x = fusion.gene.id, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 4, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="",
       fill = "Read Support") +
  theme_minimal() +
  facet_nested(Sequence_Identity ~ depth + recall_category , scales = "free_x", space = "free")+
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 50)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.spacing.y = unit(0.2, "lines"))  # reduce vertical space
ggsave("supp_Figure11_A.pdf", plot = supp_Figure11_A,  width = 400, height = 210 , unit ="mm")
ggsave("supp_Figure11_A.png", plot = supp_Figure11_A,  width = 400, height = 210 , unit ="mm")


supp_Figure11_B <-ggplot(filter(fusions_in_n_algorithms, depth == "100GB", control == "positive", overlap >= 2, recall_category == "True_Recall"), aes(x = fusion.gene.id, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 3, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="",
       fill = "Read Support") +
  theme_minimal() +
  facet_nested(Sequence_Identity ~ depth + recall_category , scales = "free_x", space = "free")+
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 50)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.spacing.y = unit(0.2, "lines"))  # reduce vertical space
ggsave("supp_Figure11_B.pdf", plot = supp_Figure11_B,  width = 400, height = 210 , unit ="mm")
ggsave("supp_Figure11_B.png", plot = supp_Figure11_B,  width = 400, height = 210 , unit ="mm")

supp_Figure11 <- supp_Figure11_A / supp_Figure11_B
ggsave("supp_Figure11.png", plot = supp_Figure11,  width = 400, height = 420 , unit ="mm")
ggsave("supp_Figure11.pdf", plot = supp_Figure11,  width = 400, height = 420 , unit ="mm")
#Supplementary figure 12####
supp_Figure12 <-  ggplot(filter(Overlap_algorithms, Sequencing_Depth == "Total", overlap > 3), 
       aes(x = fusionGeneID, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 4, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="", fill = "") +
  theme_minimal() +
  facet_grid(Library ~ Cell_Lines, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 30)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("supp_Figure12.png", plot = supp_Figure12,  width = 297, height = 210 , unit ="mm")
ggsave("supp_Figure12.pdf", plot = supp_Figure12,  width = 297, height = 210 , unit ="mm")
#Supplementary figure 13####
supp_Figure13 <-ggplot(filter(absolute_Breakpoint_Accuracy), aes(y = values_vector, x = Breakpoint_number, color = Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina()+
  theme_minimal()+ 
  scale_y_continuous(trans='log1p', labels = label_log(digits = 2),
                     breaks = c(0, 1, 10, 100, 1000, 10000, 100000, 1000000, 1000000, 10000000, 100000000, 1000000000))+
  labs(y = "Absolute distance from simulated breakpoint", x= "")+
  facet_grid(seq_depth ~ seq_id)+ 
  theme(panel.grid.minor = element_blank(),   # remove minor log lines
        panel.grid.major.x = element_blank(),  # optional: remove vertical major linesplot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("supp_Figure13.png", plot = supp_Figure13,  width = 297, height = 210 , unit ="mm")
ggsave("supp_Figure13.pdf", plot = supp_Figure13,  width = 297, height = 420 , unit ="mm")
#Supplementary figure 14 ####
supp_Figure14A <- ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%"),
                        aes(y = spanning.reads, x = classification, colour =recall_category ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL’s classification of False, Partial and True fusion calls ",
        title = "Sequencing Depth 100Gb , 95% Mean Sequence Identity",
       x = "", y ="supporting reads")
supp_Figure14B <- ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%", fusionType=="read_through"),
       aes(y = spanning.reads, x = classification, colour =fusionType ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL's classication of simulated read-throughs",
       x="",y ="supporting reads")
supp_Figure14C <-ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%", classification=="PotentialReadThrough"),
                        aes(y = spanning.reads, x = classification, colour = recall_category  ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL's 'Potential Readthrough' recall classification", x="", y="")+
  ylim(c(0,20))
supp_Figure14 <- supp_Figure14A/(supp_Figure14B + supp_Figure14C)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")),
                  tag_levels = "A")
  
ggsave("supp_Figure14.png", plot = supp_Figure14,  width = 297, height = 210 , unit ="mm")
ggsave("supp_Figure14.pdf", plot = supp_Figure14,  width = 297, height = 210 , unit ="mm")

#Supplementary figure 15 ####
ggplot(filter(Genion_Fffscores, control =="positive"))+
  geom_point(aes(x = V4, y= V3.x, colour = recall_category))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10(breaks = 10^(-5:3))+
  labs(x="FiN score", y="ffigf score")+ 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
#Supplementary figure 16 ####
ggplot(filter(Genion_Fffscores, control =="positive", fusionType=="read_through"))+
  geom_point(aes(x = V4, y= V3.x, colour = V7))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10(breaks = 10^(-5:5))+
  labs(x="FiN score", y="ffigf score")+ 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
#Supplementary figure 17 ####
supp_Figure17 <- ggplot(filter(obv_false, NumSupp >=2, Sequencing_Depth == c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb")), aes(x = Sequencing_Depth, fill = fusionType)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "False calls by each algorithm", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(Library+Cell_Lines~Algorithm) +
  scale_y_log10()+
  labs(fill = "Fusion Type") 
ggsave("supp_Figure17.png", plot = supp_Figure17,  width = 210, height =297  , unit ="mm")
ggsave("supp_Figure17.pdf", plot = supp_Figure17,  width = 210, height =297  , unit ="mm")
