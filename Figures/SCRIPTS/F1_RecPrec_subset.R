#F1, Precision & Recall Subsets
F1_subset<- ggplot(filter(readsupp_filtering, control == "positive", depth == "100GB", Sequence_Identity == "95%"))+
  geom_line(aes(x = minimum_read_support, y = F1, colour = Algorithm))+
  scale_x_log10()+
  labs(x = "Read Support",
       y = "F1")+
  ylim(c(0,1.00)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
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

Figure3A <-((Precision_subset /  Recall_subset) | F1_subset) + 
  plot_layout(guides = 'keep')   