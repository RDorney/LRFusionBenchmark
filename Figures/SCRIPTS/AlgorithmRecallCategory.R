Genion_Fffscores <- Annot_Genion_Sim[c(1:27, 37, 44:46)]%>%
  unique()%>%  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall"))%>%
  dplyr::select(Source, V4, V3.x, recall_category, control, depth, Sequence_Identity, fusionType, V7) %>% unique() #remove duplicate rows of true fusions

write_tsv(Genion_Fffscores, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/Genion_Fffscores.tsv")

JAFFAL_scoring <- Annot_JAFFAL_Sim[c(1:25, 43, 50:52)] %>% unique() %>%  
  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall")) %>%
  dplyr::select(sample, contig, fusion.genes, spanning.reads, classification, recall_category, control, depth, Sequence_Identity, fusionType) %>% unique()

write_tsv(JAFFAL_scoring, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/JAFFAL_scoring.tsv")

GFSeeker_scoring <- Annot_GFSeeker_Sim[c(1:17, 39, 46, 47)] %>% 
  unique() %>%  
  mutate(recall_category = case_when(
    grepl("false", fusionType) ~ "False_Call",
    grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
    TRUE ~ "True_Recall")) %>%
  dplyr::select(Source, 'supporting reads information', fusion.gene.id, 
                "support num", "rank class", 
                recall_category, control, depth, Sequence_Identity, fusionType) %>% unique()

write_tsv(GFSeeker_scoring, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/GFSeeker_scoring.tsv")

#Check if JAFFAL's classification correlates to recall category. 
#Make Supplementary Figure 15
Supp_Figure15A <- ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%"),
                         aes(y = spanning.reads, x = classification, colour =recall_category ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL’s classification of False, Partial and True fusion calls ",
       title = "Sequencing Depth 100Gb , 95% Mean Sequence Identity",
       x = "", y ="supporting reads")

Supp_Figure15B <- ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%", fusionType=="read_through"),
                         aes(y = spanning.reads, x = classification, colour =fusionType ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL's classication of simulated read-throughs",
       x="",y ="supporting reads")
supp_Figure15C <-ggplot(filter(JAFFAL_scoring, control =="positive", depth == "100GB", Sequence_Identity == "95%", classification=="PotentialReadThrough"),
                        aes(y = spanning.reads, x = classification, colour = recall_category  ))+
  geom_boxplot(outlier.shape = NA)+
  geom_sina()+  
  theme_minimal()+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  labs(subtitle = "JAFFAL's 'Potential Readthrough' recall classification", x="", y="")+
  ylim(c(0,20))
supp_Figure15 <- Supp_Figure15A/(Supp_Figure15B + supp_Figure15C)+
  plot_annotation(theme = theme(plot.title = element_text(size = 14, face = "bold")),
                  tag_levels = "A")
supp_Figure15
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure15.png", plot = supp_Figure15,  width = 297, height = 210 , unit ="mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure15.pdf", plot = supp_Figure15,  width = 297, height = 210 , unit ="mm")

#check if Genion's FiN or ffigf score correlate to recall category (True, Partial or False call) or other variables. 
#Make Supplementary Figure 16 & 17
#supplementary figure 16
supp_Figure16 <- ggplot(filter(Genion_Fffscores, control =="positive"))+
  geom_point(aes(x = V4, y= V3.x, colour = recall_category))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10(breaks = 10^(-5:3))+
  labs(x="FiN score", y="ffigf score")+ 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure16.png", plot = supp_Figure16,  width = 297, height = 210 , unit ="mm")
#supplementary figure 17
supp_Figure17 <- ggplot(filter(Genion_Fffscores, control =="positive", fusionType=="read_through"))+
  geom_point(aes(x = V4, y= V3.x, colour = V7))+
  facet_grid(depth ~ Sequence_Identity)+
  scale_x_log10(breaks = 10^(-5:5))+
  labs(x="FiN score", y="ffigf score")+ 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure17.png", plot = supp_Figure17,  width = 297, height = 210 , unit ="mm")
