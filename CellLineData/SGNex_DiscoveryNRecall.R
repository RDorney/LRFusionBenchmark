#####################################################
# Check Discovery & Recall of Fusions in K562 & MCF7
#####################################################

#####################
# JAFFAL
#####################
disc_JAFFAL_SGNex <- left_join(JAFFAL_SGNex_Annot, known_fusions_manual_annotation, 
                               by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_JAFFAL_SGNex$fusion.genes[disc_JAFFAL_SGNex$discovery=="known"] %>% 
  unique()

disc_JAFFAL_SGNex <-disc_JAFFAL_SGNex %>%  mutate(discovery = case_when(
  fusion.genes %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery)) %>%  
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

disc_JAFFAL_SGNex_3Gene <- JAFFAL_SGNex_Annot_3Gene %>%  
  mutate(discovery ="novel", fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y, "::", ensembl_gene_id)) %>% 
  dplyr::rename(fusion.genes = Fusion, spanning.reads = Reads, classification = Classification)

SGNex_summary_JAFFAL_df <- disc_JAFFAL_SGNex%>% 
  full_join(disc_JAFFAL_SGNex_3Gene) %>%
  group_by(across(all_of(colnames(JAFFAL_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_JAFFAL_df <- SGNex_summary_JAFFAL_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# Genion
#####################
disc_Genion_SGNex <- left_join(Genion_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("Gene1" = "ensembl_gene_id.x", 
                                     "Gene2" = "ensembl_gene_id.y", 
                                     "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_Genion_SGNe$V1[disc_Genion_SGNe$discovery=="known"] %>% unique()
disc_Genion_SGNex <- disc_Genion_SGNex%>%  mutate(discovery = case_when(
  V1 %in% check ~ "known",  # if the fusion is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_Genion_df <- disc_Genion_SGNex%>%
  group_by(across(all_of(colnames(Genion_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_Genion_df <- SGNex_summary_Genion_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# FusionSeeker
#####################
disc_FusionSeeker_SGNex <- left_join(FusionSeeker_SGNex_Annot, known_fusions_manual_annotation, 
                                    by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_FusionSeeker_SGNex$fusionGeneID[disc_FusionSeeker_SGNex$discovery=="known"] %>% 
  unique()
disc_FusionSeeker_SGNex<-disc_FusionSeeker_SGNex%>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%  
  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_FusionSeeker_df <- disc_FusionSeeker_SGNex%>%
  group_by(across(all_of(colnames(FusionSeeker_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_FusionSeeker_df <- SGNex_summary_FusionSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# LongGF
#####################
disc_LongGF_SGNex <- left_join(LongGF_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("ensembl_gene_id.x", "ensembl_gene_id.y", 
                                     "Cell_Lines"="Cell_Line")) %>% 
  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_LongGF_SGNex$fusionGeneID[disc_LongGF_SGNex$discovery=="known"] %>% 
  unique()
disc_LongGF_SGNex<-disc_LongGF_SGNex%>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_LongGF_df <- disc_LongGF_SGNex%>%
  group_by(across(all_of(colnames(LongGF_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_LongGF_df <- SGNex_summary_LongGF_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# CTAT-LR-Fusion
#####################
disc_CTATLR_SGNex <- left_join(CTATLR_SGNex_Annot, known_fusions_manual_annotation, 
                               by = c("ensembl_gene_id.x", "ensembl_gene_id.y", 
                                      "Cell_Lines"="Cell_Line")) %>% 
  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_CTATLR_SGNex$fusionGeneID[disc_CTATLR_SGNex$discovery=="known"] %>% 
  unique()
disc_CTATLR_SGNex<-disc_CTATLR_SGNex%>% 
  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_CTATLR_df <- disc_CTATLR_SGNex%>%
  group_by(across(all_of(colnames(CTATLR_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_CTATLR_df <- SGNex_summary_CTATLR_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# GFSeeker
#####################
disc_GFSeeker_SGNex<- left_join(GFSeeker_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("ensembl_gene_id.x", "ensembl_gene_id.y", 
                                     "Cell_Lines"="Cell_Line")) %>% 
  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_GFSeeker_SGNe$fusionGeneID[disc_GFSeeker_SGNe$discovery=="known"] %>% 
  unique()
disc_GFSeeker_SGNex<-disc_GFSeeker_SGNex%>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_GFSeeker_df <- disc_GFSeeker_SGNex%>%
  group_by(across(all_of(colnames(GFSeeker_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_GFSeeker_df <- SGNex_summary_GFSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

######################################################
# Combine results from different tools for comparison
######################################################
super_SGNex_summary <- rbind(SGNex_summary_FusionSeeker_df, SGNex_summary_GFSeeker_df, 
                             SGNex_summary_Genion_df, SGNex_summary_JAFFAL_df, 
                             SGNex_summary_LongGF_df, SGNex_summary_CTATLR_df)
super_SGNex_summary$Sequencing_Depth<- factor(super_SGNex_summary$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
super_SGNex_summary$Library<- factor(super_SGNex_summary$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
super_SGNex_summary %>% filter(Sequencing_Depth != "Total")

#Figure 5 ####
# Plot for K562 with smaller y-axis limit
p_k562 <- ggplot(filter(super_SGNex_summary, Sequencing_Depth != "Total", Cell_Lines == "K562")) +
  geom_point(aes(x = novelnumber, y = knownnumber, colour = Algorithm, size = Sequencing_Depth)) +
  facet_grid(Library ~ Cell_Lines)+ 
  theme(legend.position = "none",
        strip.text.y = element_blank(),         # removes row facet labels (Library)
        strip.background.y = element_blank(),    # removes grey background for row labels
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  scale_x_continuous(trans='log1p', 
                     breaks = c(0, 1, 10, 100, 
                                1000, 10000),
                     labels = c("0", "1", "10", expression(10^2),
                                expression(10^3), expression(10^4)))+
  labs(x = "Novel Fusions", y = "Known Fusions") +
  ylim(0, 6)  # Adjust this limit to your preference

# Plot for MCF7 with larger y-axis limit
p_mcf7 <- ggplot(filter(super_SGNex_summary, Sequencing_Depth != "Total", Cell_Lines == "MCF7")) +
  geom_point(aes(x = novelnumber, y = knownnumber, colour = Algorithm, size = Sequencing_Depth)) +
  facet_grid(Library ~ Cell_Lines) +
  labs(x = "Novel Fusions", y = "")+ 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  scale_x_continuous(trans='log1p', 
                     breaks = c(0, 1, 10, 100, 
                                1000, 10000),
                     labels = c("0", "1", "10", expression(10^2),
                                expression(10^3), expression(10^4)))+
  ylim(0, 55)  # Adjust as needed

#combine the two plots
Figure5 <- p_k562 + p_mcf7
ggsave("Figure5_paper.pdf", plot = Figure5, width = 210, height =200 , units = "mm")
ggsave("Figure5_paper.jpeg", plot = Figure5, width = 210, height =200 , units = "mm")

disc_JAFFAL_SGNex<- left_join(JAFFAL_SGNex_Annot, DepMap_specificCells[c(2,5,7,8,12,16)], by = c(
  "Gene1" = "Gene1_name",
  "Gene2" = "Gene2_name",
  "ensembl_gene_id.x" = "ensembl_gene_id_Gene1",
  "ensembl_gene_id.y" = "ensembl_gene_id_Gene2",
  "Cell_Lines" = "StrippedCellLineName"
)) %>%
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_JAFFAL_SGNe$fusion.genes[disc_JAFFAL_SGNe$discovery=="known"] %>% unique()
disc_JAFFAL_SGNex<-disc_JAFFAL_SGNex%>%  mutate(discovery = case_when(
  fusion.genes %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery)) %>%  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
disc_JAFFAL_SGNex_3Gene <- JAFFAL_SGNex_Annot_3Gene %>%  mutate(discovery ="novel", fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y, "::", ensembl_gene_id)) %>% rename(fusion.genes = Fusion, spanning.reads = Reads, classification = Classification)

SGNex_summary_JAFFAL_df <- disc_JAFFAL_SGNex%>% full_join(disc_JAFFAL_SGNex_3Gene) %>%
  group_by(across(all_of(colnames(JAFFAL_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_JAFFAL_df <- SGNex_summary_JAFFAL_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_Genion_SGNex<- left_join(Genion_SGNex_Annot, DepMap_specificCells[c(2,3,7,8,12,16)], by = c("Gene1" = "ensembl_gene_id.x", "Gene2" = "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_Genion_SGNe$V1[disc_Genion_SGNe$discovery=="known"] %>% unique()
disc_Genion_SGNex<-disc_Genion_SGNex%>%  mutate(discovery = case_when(
  V1 %in% check ~ "known",  # if the fusion is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_Genion_df <- disc_Genion_SGNex%>%
  group_by(across(all_of(colnames(Genion_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_Genion_df <- SGNex_summary_Genion_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_FusionSeeker_SGNex<- left_join(FusionSeeker_SGNex_Annot, DepMap_specificCells[c(2,3,7,8,12,16)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_FusionSeeker_SGNe$fusionGeneID[disc_FusionSeeker_SGNe$discovery=="known"] %>% unique()
disc_FusionSeeker_SGNex<-disc_FusionSeeker_SGNex%>% mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_FusionSeeker_df <- disc_FusionSeeker_SGNex%>%
  group_by(across(all_of(colnames(FusionSeeker_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_FusionSeeker_df <- SGNex_summary_FusionSeeker_df %>% group_by(Source, Cell_Lines, 
                                                                                  Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_LongGF_SGNex<- left_join(LongGF_SGNex_Annot, DepMap_specificCells[c(2,3,7,8,12,16)], 
                              by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_LongGF_SGNe$fusionGeneID[disc_LongGF_SGNe$discovery=="known"] %>% 
  unique()
disc_LongGF_SGNex<-disc_LongGF_SGNex%>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_LongGF_df <- disc_LongGF_SGNex%>%
  group_by(across(all_of(colnames(LongGF_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
SGNex_summary_LongGF_df <- SGNex_summary_LongGF_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

super_SGNex_summary <- rbind(SGNex_summary_FusionSeeker_df, SGNex_summary_Genion_df, SGNex_summary_JAFFAL_df, SGNex_summary_LongGF_df)
# Convert Sequencing_Depth into an ordered factor
super_SGNex_summary$Sequencing_Depth<- factor(super_SGNex_summary$Sequencing_Depth, 
                                              levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
# Convert Library type into an ordered factor
super_SGNex_summary$Library<- factor(super_SGNex_summary$Library, 
                                     levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
# Quick check: filter out "Total" depth entries
super_SGNex_summary %>% filter(Sequencing_Depth != "Total")

SGNex_novelvsknown <- ggplot(filter(super_SGNex_summary, Sequencing_Depth != "Total"))+
  geom_point(aes(x = novelnumber , y =knownnumber, colour = Algorithm, shape = Sequencing_Depth, ))+
  facet_grid(Library~Cell_Lines, scales = "free_y")+
  labs(x="Novel Fusions", y = "Known Fusions") 
ggsave("SGNex_novelvsknown.pdf", plot = SGNex_novelvsknown, width = 297, height = 210, units = "mm")