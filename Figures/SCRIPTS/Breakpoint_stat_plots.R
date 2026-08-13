absolute_Breakpoint_Accuracy$Distance_from_Simulated_Breakpoint <-
  as.numeric(as.character(absolute_Breakpoint_Accuracy$Distance_from_Simulated_Breakpoint))

#Statistical test: Check if data is normally distributed
signed_log1p <- function(x) sign(x) * log1p(abs(x))
breakpoint_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%")
breakpoint_stat_data$Signed_log_values <- signed_log1p(breakpoint_stat_data$Distance_from_Simulated_Breakpoint)
by(breakpoint_stat_data$Signed_log_values,
   breakpoint_stat_data$Algorithm,
   shapiro.test)
breakpoint1_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%", Breakpoint_number == "Breakpoint_1")
breakpoint1_stat_data$Signed_log_values <- signed_log1p(breakpoint1_stat_data$Distance_from_Simulated_Breakpoint)
by(breakpoint1_stat_data$Signed_log_values,
   breakpoint1_stat_data$Algorithm,
   shapiro.test)
breakpoint2_stat_data <- filter(absolute_Breakpoint_Accuracy, seq_depth == "100GB", seq_id == "95%", Breakpoint_number == "Breakpoint_2")
breakpoint2_stat_data$Signed_log_values <- signed_log1p(breakpoint2_stat_data$Distance_from_Simulated_Breakpoint)
by(breakpoint2_stat_data$Signed_log_values,
   breakpoint2_stat_data$Algorithm,
   shapiro.test)

#Data is not normally distributed, so we will do a Kruskal Wallis Test followed by a pairwise test.
kruskal.test(Signed_log_values ~ Algorithm, data = breakpoint_stat_data)
kruskal.test(Signed_log_values ~ Algorithm, data = breakpoint1_stat_data)
kruskal.test(Signed_log_values ~ Algorithm, data = breakpoint2_stat_data)

library(dplyr)
library(rstatix)

one_sample_tests <- absolute_Breakpoint_Accuracy %>%
  filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Breakpoint_number, Algorithm) %>%
  wilcox_test(Distance_from_Simulated_Breakpoint ~ 1, mu = 0) %>%
  mutate(p.signif = p_format(p, digits = 2))%>%
  mutate(y.position = 7.5)

##########################################################
# check if breakpoints are significantly different from 0
##########################################################
wilcox_table<- absolute_Breakpoint_Accuracy %>%
  dplyr::filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Breakpoint_number, Algorithm) %>%
  dplyr::reframe(
    # Perform the Wilcoxon Signed-Rank Test against 0
    test = list(wilcox.test(Distance_from_Simulated_Breakpoint, 
                            mu = 0, 
                            alternative = "two.sided", 
                            exact = FALSE)),
    # Use broom::tidy to turn the test object into a row of data
    tidy(test[[1]]) 
  ) %>%
  # Clean up the output table
  dplyr::select(Breakpoint_number, Algorithm, statistic, p.value, method) %>%
  mutate(significance = dplyr::case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

library(infer)
bootstrap_table <-absolute_Breakpoint_Accuracy %>%
  dplyr::filter(seq_depth == "100GB", seq_id == "95%") %>%
  group_by(Tool_Breakpoint, Algorithm) %>%
  # We use group_modify to run the bootstrap for EACH group
  dplyr::group_modify(~ {
    .x %>%
      specify(response = Distance_from_Simulated_Breakpoint) %>%
      generate(reps = 10000, type = "bootstrap") %>%
      calculate(stat = "median") %>%
      get_confidence_interval(level = 0.95, type = "percentile")
  })
#CTATLRB2, FSB2, GFSB2, GB2 (very precise), LB2 span 0, therefore n.s

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
  labs(y = "Distance from simulated breakpoint", x= "")+
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
