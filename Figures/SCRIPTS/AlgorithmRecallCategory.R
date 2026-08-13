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

