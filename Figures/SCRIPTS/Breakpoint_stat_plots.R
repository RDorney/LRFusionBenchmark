#############################
# Plot Breakpoint distances #
#############################
# Make Figure 4 &
# Supplementary Figure 13
absolute_Breakpoint_Accuracy$Distance_from_Simulated_Breakpoint <-
  as.numeric(as.character(absolute_Breakpoint_Accuracy$Distance_from_Simulated_Breakpoint))

#Statistical test: Check if data is normally distributed
breakpoint_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%")
breakpoint_stat_data$Log1p_values <- log1p(breakpoint_stat_data$values_vector)
by(breakpoint_stat_data$Log1p_values,
   breakpoint_stat_data$Algorithm,
   shapiro.test)
breakpoint1_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%", Breakpoint_number == "Breakpoint_1")
breakpoint1_stat_data$Log1p_values <- log1p(breakpoint1_stat_data$values_vector)
by(breakpoint1_stat_data$Log1p_values,
   breakpoint1_stat_data$Algorithm,
   shapiro.test)
breakpoint2_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%", Breakpoint_number == "Breakpoint_2")
breakpoint2_stat_data$Log1p_values <- log1p(breakpoint2_stat_data$values_vector)
by(breakpoint2_stat_data$Log1p_values,
   breakpoint2_stat_data$Algorithm,
   shapiro.test)

ggplot(breakpoint_stat_data,
       aes(sample = Log1p_values)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_grid(Breakpoint_number ~ Algorithm) +
  theme_minimal() +
  labs(title = "Q–Q Plots by Algorithm and Breakpoint",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
#Data is not normally distributed, so we will do a Kruskal Wallis Test followed by a pairwise test.
kruskal.test(Log1p_values ~ Algorithm, data = breakpoint_stat_data)
kruskal.test(Log1p_values ~ Algorithm, data = breakpoint1_stat_data)
kruskal.test(Log1p_values ~ Algorithm, data = breakpoint2_stat_data)

#Figure 4 ####
Figure4 <- ggplot(filter(absolute_Breakpoint_Accuracy, 
                         seq_depth == "100GB", seq_id == "95%"), 
                  aes(y = Distance_from_Simulated_Breakpoint, x = Breakpoint_number, 
                      color = Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina(maxwidth = 0.8)+
  theme_minimal()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid.minor = element_blank(),   # remove minor log lines
        panel.grid.major.x = element_blank(),  # optional: remove vertical major lines
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))+ 
  scale_y_continuous(transform = scales::pseudo_log_trans(base=10) , 
                     breaks = c(-1000000000, -100000000, 
                                -10000000, -1000000, -100000, 
                                -10000, -1000, -100, 
                                -10, -1, 
                                0, 
                                1, 10, 
                                100, 1000, 10000, 
                                100000, 1000000, 10000000, 
                                100000000, 1000000000),
                     labels = c(expression(-10^9), expression(-10^8), 
                                expression(-10^7), expression(-10^6), expression(-10^5), 
                                expression(-10^4), expression(-10^3), expression(-10^2),
                                expression(-10^1), "-1", 
                                "0", 
                                "1", expression(10^1), 
                                expression(10^2), expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6), expression(10^7), 
                                expression(10^8), expression(10^9))
  )+   
  scale_x_discrete(labels = c("Breakpoint 1", "Breakpoint 2"))+
  labs(y = "Distance from simulated breakpoint",
       x= "")+
  stat_compare_means(method = "kruskal.test", aes(group = Algorithm), label.y = 8.5 , show.legend = FALSE)+
  geom_pwc(method = "dunn_test", p.adjust.method = "bonferroni", hide.ns = TRUE, label = "p.adj.signif",
           y.position = c(6.5, 7.5 , 5.5, 6.5 ,7.5 ), show.legend = FALSE)

ggsave("~/LongReadFusionCallerBenchmark/Figures/Figure4_paper.pdf", plot = Figure4, width = 297, height = 210, units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/Figure4_paper.png", plot = Figure4)

supp_Figure13 <-ggplot(filter(absolute_Breakpoint_Accuracy), aes(y = Distance_from_Simulated_Breakpoint, x = Breakpoint_number, color = Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_sina()+
  theme_minimal()+ 
  scale_y_continuous(transform = scales::pseudo_log_trans(base=10) , 
                     breaks = c(-1000000000, -100000000, 
                                -10000000, -1000000, -100000, 
                                -10000, -1000, -100, 
                                -10, -1, 
                                0, 
                                1, 10, 
                                100, 1000, 10000, 
                                100000, 1000000, 10000000, 
                                100000000, 1000000000),
                     labels = c(expression(-10^9), expression(-10^8), expression(-10^7), expression(-10^6), expression(-10^5), expression(-10^4), expression(-10^3), expression(-10^2),
                                expression(-10^1), "-1", 
                                "0", 
                                "1", expression(10^1), 
                                expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8), expression(10^9))
  )+   
  scale_x_discrete(labels = c("Breakpoint 1", "Breakpoint 2", "Breakpoint 3", "Breakpoint 4"))+
  labs(y = "Absolute distance from simulated breakpoint", x= "")+
  stat_compare_means(method = "kruskal.test", aes(group = Algorithm), label.y = 8.5 , show.legend = FALSE, label = "p.format",)+
  geom_pwc(method = "dunn_test", p.adjust.method = "bonferroni", hide.ns = TRUE, label = "p.adj.signif",
           y.position = c(6.5, 7.5 , 5.5, 6.5 ,7.5 ), show.legend = FALSE)+
  facet_grid(seq_depth ~ seq_id)+ 
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid.minor = element_blank(),   # remove minor log lines
        panel.grid.major.x = element_blank(),  # optional: remove vertical major linesplot.title = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size = 9.5),
        axis.title = element_text(size = 12))
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure13.png", plot = supp_Figure13,  width = 300, height = 420 , unit ="mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure13.pdf", plot = supp_Figure13,  width = 297, height = 420 , unit ="mm")