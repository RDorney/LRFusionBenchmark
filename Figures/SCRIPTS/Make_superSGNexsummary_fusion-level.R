
#super_SGNex_summary_fusion_level <-
  
  disc_JAFFAL_SGNex%>% 
  full_join(disc_JAFFAL_SGNex_3Gene) %>%
  dplyr::group_by(across(c("sample","fusion.genes","chrom1","chrom2",
                           "Source","Cell_Lines","Algorithm",
                           "Sequencing_Depth","Library","fusionType",
                           "discovery","fusionGeneID" ))) %>%
  dplyr::summarise(read_supp = sum(spanning.reads), .groups = "drop")%>%
  mutate(
    knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
    novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0)) %>%
  group_by(fusionGeneID, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, fusionGeneID)%>% 
  summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))%>%View()

  disc_Genion_SGNex%>%
    dplyr::group_by(across(c("V1","chr1","chr2",
                             "Source","Cell_Lines","Algorithm",
                             "Sequencing_Depth","Library","fusionType",
                             "discovery","fusion.gene.id"))) %>%
    dplyr::summarise(read_supp = sum(V5), .groups = "drop") %>%
    mutate(
      knowncount = case_when(discovery == "known" ~ 1, TRUE ~ 0),
      novelcount = case_when(discovery == "novel" ~ 1, TRUE ~ 0))%>% 
    group_by(Source, Cell_Lines, Sequencing_Depth, Algorithm, Library, fusionType, fusionGeneID) %>%
    summarise(knownnumber = sum(knowncount), novelnumber= sum(novelcount))
  
  disc_FusionSeeker_SGNex%>%
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
  
  disc_LongGF_SGNex%>%
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
  
  disc_CTATLR_SGNex%>%
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
  
  disc_GFSeeker_SGNex%>%
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