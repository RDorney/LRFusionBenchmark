
JAF_SGNex_summary_fusion_level <- disc_JAFFAL_SGNex%>% 
  full_join(disc_JAFFAL_SGNex_3Gene) %>%
  dplyr::group_by(across(c("sample","fusion.genes","chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID" ))) %>%
  dplyr::summarise(read_supp = sum(spanning.reads), .groups = "drop")%>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0)) %>%
  group_by(fusionGeneID, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType)%>% 
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount)) 

Gen_SGNex_summary_fusion_level <- disc_Genion_SGNex%>%
    dplyr::group_by(across(c("V1","chr1","chr2",
                             "Source","Cell_Lines","Algorithm",
                             "Sequencing_Depth","Library","fusionType",
                             "discovery","fusion.gene.id"))) %>%
    dplyr::summarise(read_supp = sum(V5), .groups = "drop") %>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
    group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, V1) %>%
    summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))%>%
    rename(fusionGeneID = V1)
  
FusionSeeker_SGNex_summary_fusion_level <- disc_FusionSeeker_SGNex%>%
    group_by(across(c("fusionGene.x","Chrom1","Chrom2",
                      "Source","Cell_Lines","Algorithm",
                      "Sequencing_Depth","Library","fusionType",
                      "discovery","fusionGeneID"))) %>%
    dplyr::summarise(read_supp = sum(NumSupp), .groups = "drop") %>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
    group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, fusionGeneID) %>%
    summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
  
LongGF_SGNex_summary_fusion_level <- disc_LongGF_SGNex%>%
    dplyr::group_by(across(c("V2","chromosome1","chromosome2",
                             "Source","Cell_Lines","Algorithm",
                             "Sequencing_Depth","Library","fusionType",
                             "discovery","fusionGeneID"))) %>%
    dplyr::summarise(read_supp = sum(V3), .groups = "drop") %>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
    group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, fusionGeneID) %>%
    summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
  
CTATLR_SGNex_summary_fusion_level <- disc_CTATLR_SGNex%>%
    dplyr::group_by(across(c("#FusionName","LeftGene","RightGene",
                             "chrom1","chrom2",
                             "Source","Cell_Lines","Algorithm",
                             "Sequencing_Depth","Library","fusionType",
                             "discovery","fusionGeneID"))) %>%
    dplyr::summarise(read_supp = sum(num_LR), .groups = "drop")%>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
    group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, fusionGeneID) %>%
    summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
  
GFSeeker_SGNex_summary_fusion_level <- disc_GFSeeker_SGNex%>%
    dplyr::group_by(across(c("gene1_name","gene2_name","chrom1","chrom2",
                             "Source","Cell_Lines","Algorithm",
                             "Sequencing_Depth","Library","fusionType",
                             "discovery","fusionGeneID"))) %>%
    dplyr::summarise(read_supp = sum(`support num`), .groups = "drop") %>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
  group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType,fusionGeneID) %>%
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))

super_SGNex_summary_fusionlevel <- rbind(JAF_SGNex_summary_fusion_level, Gen_SGNex_summary_fusion_level, 
                                         FusionSeeker_SGNex_summary_fusion_level, LongGF_SGNex_summary_fusion_level, 
                                         CTATLR_SGNex_summary_fusion_level, GFSeeker_SGNex_summary_fusion_level)

write_tsv(super_SGNex_summary_fusionlevel, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/fusionlevel_SGNex_summary.tsv.gz")

library(dplyr)

JAF_SGNex_summary_library_level <- disc_JAFFAL_SGNex %>% 
  full_join(disc_JAFFAL_SGNex_3Gene) %>%
  dplyr::group_by(sample, fusion.genes, chrom1, chrom2,
                  Source, Cell_Lines, Algorithm,
                  Sequencing_Depth, Library, fusionType,
                  discovery, fusionGeneID) %>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber) %>%
  unique()
# View(JAF_SGNex_summary_library_level)

Gen_SGNex_summary_library_level <- disc_Genion_SGNex%>%
  dplyr::group_by(across(c("V1","chr1","chr2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusion.gene.id"))) %>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  rename(fusionGeneID = V1)%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber)%>%
  unique()

FusionSeeker_SGNex_summary_library_level <- disc_FusionSeeker_SGNex%>%
  group_by(across(c("fusionGene.x","Chrom1","Chrom2",
                    "Source","Cell_Lines","Algorithm",
                    "Sequencing_Depth","Library","fusionType",
                    "discovery","fusionGeneID"))) %>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber)%>%
  unique()

LongGF_SGNex_summary_library_level <- disc_LongGF_SGNex%>%
  dplyr::group_by(across(c("V2","chromosome1","chromosome2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID")))%>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber)%>%
  unique()

CTATLR_SGNex_summary_library_level <- disc_CTATLR_SGNex%>%
  dplyr::group_by(across(c("#FusionName","LeftGene","RightGene",
                           "chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID"))) %>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber)%>%
  unique()

GFSeeker_SGNex_summary_library_level <- disc_GFSeeker_SGNex%>%
  dplyr::group_by(across(c("gene1_name","gene2_name","chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID"))) %>%
  dplyr::summarise(
    knownnumber = sum(discovery == "known"), 
    novelnumber = sum(discovery == "novel"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    knownnumber = case_when(knownnumber != 0 ~ 1, TRUE ~ 0),
    novelnumber = case_when(novelnumber != 0 ~ 1, TRUE ~ 0)
  )%>%
  select(Cell_Lines, Algorithm, fusionType, fusionGeneID, knownnumber, novelnumber)%>%
  unique()

super_SGNex_summary_library_level <- rbind(JAF_SGNex_summary_library_level, Gen_SGNex_summary_library_level, 
                                         FusionSeeker_SGNex_summary_library_level, LongGF_SGNex_summary_library_level, 
                                         CTATLR_SGNex_summary_library_level, GFSeeker_SGNex_summary_library_level)

# View(super_SGNex_summary_library_level)
write_tsv(super_SGNex_summary_library_level, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/library_level_SGNex_summary.tsv.gz")
