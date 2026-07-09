#F1, Precision & Recall Subsets
F1_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "F1")+
  ylim(c(0,0.40)) +
  scale_colour_npg()+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  theme_bw()
Precision_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = Precision, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "", y = "Precision") +
  ylim(c(0,1.0))+
  scale_colour_npg()+
  theme_bw()+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) 
Recall_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = Recall, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "Recall")+
  ylim(c(0,0.40))+
  scale_colour_npg()+
  theme_bw()+
  theme(plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

Figure3A <-(F1_subset / Precision_subset /  Recall_subset) + 
  plot_layout(guides = 'keep')   
