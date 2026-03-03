#Supplementary figure 8####
Recall_Precision <- ggplot(filter(stat_summary, control == 'positive'))+
  geom_point(aes(x = Precision, y= Recall, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_colour_npg()+
  theme_bw()+
  labs(x = "Precision", y = "Recall")+
  ylim(c(0,.5))+
  xlim(c(0,.6)) +
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/Recall-Precision.pdf", plot = Recall_Precision, width = 297, height = 210, units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/Recall-Precision.png", plot = Recall_Precision, width = 297, height = 210, units = "mm")
#Supplementary figure 9####
ReadSupport_Precision <- ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "Precision") +
  scale_colour_npg()+
  theme_bw()+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/ReadSupport-Precision.pdf", plot = ReadSupport_Precision, width = 297, height = 210, units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/ReadSupport-Precision.png", plot = ReadSupport_Precision, width = 297, height = 210, units = "mm")
#Supplementary figure 10 ####
ReadSupport_F1Score <- ggplot(filter(readsupp_filtering, control == "positive"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "F1 Score")+
  ylim(c(0,0.40))+
  scale_colour_npg()+
  theme_bw()+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/ReadSupport-F1Score.pdf", plot = ReadSupport_F1Score, width = 297, height = 210, units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/FigureFiles/ReadSupport-F1Score.png", plot = ReadSupport_F1Score, width = 297, height = 210, units = "mm")
