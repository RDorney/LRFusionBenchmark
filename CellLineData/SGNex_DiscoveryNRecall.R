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
  dplyr::group_by(across(c("sample","fusion.genes","chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID" ))) %>%
  dplyr::summarise(read_supp = sum(spanning.reads), .groups = "drop")%>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))

total_counts_SGNex_JAFFAL_df <- SGNex_summary_JAFFAL_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_JAFFAL_df <- SGNex_summary_JAFFAL_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# Genion
#####################
disc_Genion_SGNex <- left_join(Genion_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("V1.1" = "ensembl_gene_id.x", 
                                     "V1.2" = "ensembl_gene_id.y", 
                                     "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_Genion_SGNex$V1[disc_Genion_SGNex$discovery=="known"] %>% 
  unique()
disc_Genion_SGNex <- disc_Genion_SGNex%>%  mutate(discovery = case_when(
  V1 %in% check ~ "known",  # if the fusion is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_Genion_df <- disc_Genion_SGNex%>%
  dplyr::group_by(across(c("V1","chr1","chr2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusion.gene.id"))) %>%
  dplyr::summarise(read_supp = sum(V5), .groups = "drop") %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
total_counts_SGNex_Genion_df <- SGNex_summary_Genion_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_Genion_df <- SGNex_summary_Genion_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# FusionSeeker
#####################
disc_FusionSeeker_SGNex <- left_join(FusionSeeker_SGNex_Annot, known_fusions_manual_annotation, 
                                     by = c("Gene1" = "ensembl_gene_id.x", 
                                            "Gene2" = "ensembl_gene_id.y", 
                                           "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))
check <- disc_FusionSeeker_SGNex$fusionGeneID[disc_FusionSeeker_SGNex$discovery=="known"] %>% 
  unique()
disc_FusionSeeker_SGNex<-disc_FusionSeeker_SGNex%>% 
  mutate(fusionGeneID = paste0(Gene1, "::", Gene2)) %>%  
  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_FusionSeeker_df <- disc_FusionSeeker_SGNex%>%
  group_by(across(c("fusionGene.x","Chrom1","Chrom2",
                  "Source","Cell_Lines","Algorithm",
                  "Sequencing_Depth","Library","fusionType",
                  "discovery","fusionGeneID"))) %>%
  dplyr::summarise(read_supp = sum(NumSupp), .groups = "drop") %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
total_counts_SGNex_FusionSeeker_df <- SGNex_summary_FusionSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_FusionSeeker_df <- SGNex_summary_FusionSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# LongGF
#####################
disc_LongGF_SGNex <- left_join(LongGF_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("ensembl_gene_id.x", 
                                     "ensembl_gene_id.y", 
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
  dplyr::group_by(across(c("V2","chromosome1","chromosome2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID"))) %>%
  dplyr::summarise(read_supp = sum(V3), .groups = "drop") %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
total_counts_SGNex_LongGF_df <- SGNex_summary_LongGF_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_LongGF_df <- SGNex_summary_LongGF_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# CTAT-LR-Fusion
#####################
disc_CTATLR_SGNex <- left_join(CTATLR_SGNex_Annot, known_fusions_manual_annotation, 
                               by = c("ensembl_gene_id.x", 
                                      "ensembl_gene_id.y", 
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
  dplyr::group_by(across(c("#FusionName","LeftGene","RightGene",
                           "chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID"))) %>%
  dplyr::summarise(read_supp = sum(num_LR), .groups = "drop")%>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
total_counts_SGNex_CTATLR_df <- SGNex_summary_CTATLR_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_CTATLR_df <- SGNex_summary_CTATLR_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

#####################
# GFSeeker
#####################
disc_GFSeeker_SGNex<- left_join(GFSeeker_SGNex_Annot, known_fusions_manual_annotation, 
                              by = c("ensembl_gene_id.x", "ensembl_gene_id.y", 
                                     "Cell_Lines"="Cell_Line")) %>% 
  mutate(discovery = if_else(is.na(discovery), "novel", discovery)) %>% 
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))
check <- disc_GFSeeker_SGNex$fusionGeneID[disc_GFSeeker_SGNex$discovery=="known"] %>% 
  unique()
disc_GFSeeker_SGNex<-disc_GFSeeker_SGNex%>%  mutate(discovery = case_when(
  fusionGeneID %in% check ~ "known",  # if fusion.genes is in the known list
  is.na(discovery) ~ "novel",         # if discovery is NA, label as novel
  TRUE ~ discovery))
SGNex_summary_GFSeeker_df <- disc_GFSeeker_SGNex%>%
  dplyr::group_by(across(c("gene1_name","gene2_name","chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID"))) %>%
  dplyr::summarise(read_supp = sum(`support num`), .groups = "drop") %>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))
total_counts_SGNex_GFSeeker_df <- SGNex_summary_GFSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
SGNex_summary_GFSeeker_df <- SGNex_summary_GFSeeker_df %>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

######################################################
# Combine results from different tools for comparison
######################################################
counts_SGNex_summary <- rbind(total_counts_SGNex_FusionSeeker_df, total_counts_SGNex_GFSeeker_df, 
                             total_counts_SGNex_Genion_df, total_counts_SGNex_JAFFAL_df, 
                             total_counts_SGNex_LongGF_df, total_counts_SGNex_CTATLR_df)
counts_SGNex_summary$Sequencing_Depth<- factor(counts_SGNex_summary$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
counts_SGNex_summary$Library<- factor(counts_SGNex_summary$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))

super_SGNex_summary <- rbind(SGNex_summary_FusionSeeker_df, SGNex_summary_GFSeeker_df, 
                             SGNex_summary_Genion_df, SGNex_summary_JAFFAL_df, 
                             SGNex_summary_LongGF_df, SGNex_summary_CTATLR_df)
super_SGNex_summary$Sequencing_Depth<- factor(super_SGNex_summary$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
super_SGNex_summary$Library<- factor(super_SGNex_summary$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
# Quick check: filter out "Total" depth entries
super_SGNex_summary %>% filter(Sequencing_Depth != "Total")
write_tsv(counts_SGNex_summary, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/counts_SGNex_summary.tsv.gz")
write_tsv(super_SGNex_summary, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/super_SGNex_summary.tsv.gz")

View(super_SGNex_summary)
