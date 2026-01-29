known_fusions_biomart$discovery <- "known"
fsep_JAFFALcl <- filter(fsep_JAFFALcl, spanning.reads >=2)
fsep_Gencl <- filter(fsep_Gencl, V5 >=2)
FusionSeeker_SGNex_Annot  <- filter(FusionSeeker_SGNex_Annot, NumSupp >= 2)
LongGF_SGNex_Annot <- filter(LongGF_SGNex_Annot, V3 >=2)

disc_JAFFALcl <- left_join(fsep_JAFFALcl, known_fusions_biomart[c(1,6,7,10:12)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_JAFFALcl$fusion.genes[disc_JAFFALcl$discovery=="known"] %>% unique()
disc_JAFFALcl <-disc_JAFFALcl %>%  mutate(discovery = case_when(
  fusion.genes %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery)) %>%  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
disc_JAFFALcl_3Gene <- fsep_JAFFALcl_3Gene %>%  mutate(discovery ="novel", fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y, "::", ensembl_gene_id)) %>% rename(fusion.genes = Fusion, spanning.reads = Reads, classification = Classification)

CellData_summary_JAFFAL_df <- disc_JAFFALcl %>% full_join(disc_JAFFALcl_3Gene) %>%
  group_by(across(all_of(colnames(JAFFAL_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_JAFFAL_df <- CellData_summary_JAFFAL_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_Genioncl <- left_join(fsep_Gencl, known_fusions_biomart[c(1,6,7,10:12)], by = c("Gene1" = "ensembl_gene_id.x", "Gene2" = "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_Genioncl$V1[disc_Genioncl$discovery=="known"] %>% unique()
disc_Genioncl <-disc_Genioncl %>%  mutate(discovery = case_when(
  V1 %in% check ~ "known",  # if the fusion is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_Genion_df <- disc_Genioncl %>%
  group_by(across(all_of(colnames(Genion_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_Genion_df <- CellData_summary_Genion_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_FusionSeekercl <- left_join(FusionSeeker_SGNex_Annot, known_fusions_biomart[c(1,6,7,10:12)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_FusionSeekercl$fusionGeneID[disc_FusionSeekercl$discovery=="known"] %>% unique()
disc_FusionSeekercl <-disc_FusionSeekercl %>% mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_FusionSeeker_df <- disc_FusionSeekercl %>%
  group_by(across(all_of(colnames(FusionSeeker_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_FusionSeeker_df <- CellData_summary_FusionSeeker_df %>% group_by(Source, Cell_Lines, 
                                                                                  Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_LongGFcl <- left_join(LongGF_SGNex_Annot, known_fusions_biomart[c(1,6,7,10:12)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_LongGFcl$fusionGeneID[disc_LongGFcl$discovery=="known"] %>% unique()
disc_LongGFcl <-disc_LongGFcl %>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_LongGF_df <- disc_LongGFcl %>%
  group_by(across(all_of(colnames(LongGF_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_LongGF_df <- CellData_summary_LongGF_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

super_CellData_summary <- rbind(CellData_summary_FusionSeeker_df, CellData_summary_Genion_df, CellData_summary_JAFFAL_df, CellData_summary_LongGF_df)
super_CellData_summary$Sequencing_Depth<- factor(super_CellData_summary$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
super_CellData_summary$Library<- factor(super_CellData_summary$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
super_CellData_summary %>% filter(Sequencing_Depth != "Total")

#Figure 5 ####
# Plot for K562 with smaller y-axis limit
p_k562 <- ggplot(filter(super_CellData_summary, Sequencing_Depth != "Total", Cell_Lines == "K562")) +
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
p_mcf7 <- ggplot(filter(super_CellData_summary, Sequencing_Depth != "Total", Cell_Lines == "MCF7")) +
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

disc_JAFFALcl <- left_join(fsep_JAFFALcl, DepMap_specificCells[c(2,5,7,8,12,16)], by = c(
  "Gene1" = "Gene1_name",
  "Gene2" = "Gene2_name",
  "ensembl_gene_id.x" = "ensembl_gene_id_Gene1",
  "ensembl_gene_id.y" = "ensembl_gene_id_Gene2",
  "Cell_Lines" = "StrippedCellLineName"
)) %>%
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_JAFFALcl$fusion.genes[disc_JAFFALcl$discovery=="known"] %>% unique()
disc_JAFFALcl <-disc_JAFFALcl %>%  mutate(discovery = case_when(
  fusion.genes %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery)) %>%  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
disc_JAFFALcl_3Gene <- fsep_JAFFALcl_3Gene %>%  mutate(discovery ="novel", fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y, "::", ensembl_gene_id)) %>% rename(fusion.genes = Fusion, spanning.reads = Reads, classification = Classification)

CellData_summary_JAFFAL_df <- disc_JAFFALcl %>% full_join(disc_JAFFALcl_3Gene) %>%
  group_by(across(all_of(colnames(JAFFAL_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_JAFFAL_df <- CellData_summary_JAFFAL_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_Genioncl <- left_join(fsep_Gencl, DepMap_specificCells[c(2,3,7,8,12,16)], by = c("Gene1" = "ensembl_gene_id.x", "Gene2" = "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_Genioncl$V1[disc_Genioncl$discovery=="known"] %>% unique()
disc_Genioncl <-disc_Genioncl %>%  mutate(discovery = case_when(
  V1 %in% check ~ "known",  # if the fusion is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_Genion_df <- disc_Genioncl %>%
  group_by(across(all_of(colnames(Genion_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_Genion_df <- CellData_summary_Genion_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_FusionSeekercl <- left_join(FusionSeeker_SGNex_Annot, DepMap_specificCells[c(2,3,7,8,12,16)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_FusionSeekercl$fusionGeneID[disc_FusionSeekercl$discovery=="known"] %>% unique()
disc_FusionSeekercl <-disc_FusionSeekercl %>% mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_FusionSeeker_df <- disc_FusionSeekercl %>%
  group_by(across(all_of(colnames(FusionSeeker_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_FusionSeeker_df <- CellData_summary_FusionSeeker_df %>% group_by(Source, Cell_Lines, 
                                                                                  Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

disc_LongGFcl <- left_join(LongGF_SGNex_Annot, DepMap_specificCells[c(2,3,7,8,12,16)], by = c("ensembl_gene_id.x", "ensembl_gene_id.y", "Cell_Lines"="Cell_Line")) %>%  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_LongGFcl$fusionGeneID[disc_LongGFcl$discovery=="known"] %>% unique()
disc_LongGFcl <-disc_LongGFcl %>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
CellData_summary_LongGF_df <- disc_LongGFcl %>%
  group_by(across(all_of(colnames(LongGF_SGNex)))) %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
CellData_summary_LongGF_df <- CellData_summary_LongGF_df %>% group_by(Source, Cell_Lines, 
                                                                      Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

super_CellData_summary <- rbind(CellData_summary_FusionSeeker_df, CellData_summary_Genion_df, CellData_summary_JAFFAL_df, CellData_summary_LongGF_df)
# Convert Sequencing_Depth into an ordered factor
super_CellData_summary$Sequencing_Depth<- factor(super_CellData_summary$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
# Convert Library type into an ordered factor
super_CellData_summary$Library<- factor(super_CellData_summary$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
# Quick check: filter out "Total" depth entries
super_CellData_summary %>% filter(Sequencing_Depth != "Total")

Cell_Line_novelvsknown <- ggplot(filter(super_CellData_summary, Sequencing_Depth != "Total"))+
  geom_point(aes(x = novelnumber , y =knownnumber, colour = Algorithm, shape = Sequencing_Depth, ))+
  facet_grid(Library~Cell_Lines, scales = "free_y")+
  labs(x="Novel Fusions", y = "Known Fusions") 
ggsave("Cell_Line_novelvsknown.pdf", plot = Cell_Line_novelvsknown, width = 297, height = 210, units = "mm")