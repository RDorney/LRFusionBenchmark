######################
# Coloured Table Grid#
######################
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

Supplementary_Figure2 <- (false_call_plot / partial_call_plot / true_call_plot)+ 
  plot_layout(heights = c(2, 1, 1)) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")))
