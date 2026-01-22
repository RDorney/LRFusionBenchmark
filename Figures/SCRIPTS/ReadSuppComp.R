###########################
# Read Support Comparison #
###########################
#comparing read support for the same fusions as called by different algorithms

Figure3B <-ggplot(filter(fusions_in_all_algorithms, depth == "100GB", Sequence_Identity == "95%", control == "positive", overlap == num_ALG), aes(x = fusion.gene.id, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 5, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="",
       fill = "Read Support") +
  theme_minimal() +
  facet_nested(~ recall_category , scales = "free_x", space = "free",
               labeller = as_labeller(c("Partial_Recall" = "Partial Recall", "True_Recall" = "True Recall")))+
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 25)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(strip.text.x = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.1, "cm"),   # reduce horizontal gap
        legend.spacing.y = unit(0.1, "cm"),   # reduce vertical gap inside legend
        legend.margin    = margin(0, 0, 0, 0),# trim around the legend box
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 10),
        axis.title = element_text(size = 12))+ plot_layout(guides = 'keep')

ggsave("~/LongReadFusionCallerBenchmark/Figures/Figure3B_paper.pdf", plot = Figure3B, width = 297, height = 210, units = "mm") 