#Consolidate data frames

dLGF<- Annot_LongGF_Sim[c(1:24, 38, 45:47)] %>% 
  unique()%>% #remove duplicate values
  dplyr::select(c(5, 13, 18, 25:28))%>%#select only columns needed for downstream analysis
  dplyr::rename(read_supp = V3.x)

dJF <- Annot_JAFFAL_Sim[c(1:25, 43, 50:52)] %>% 
  unique()%>% #remove duplicate values
  dplyr::select(c(13, 22, 23, 26:29)) %>% #select only columns needed for downstream analysis
  dplyr::rename(read_supp = spanning.reads)

dJF3 <- Annot_JAFFAL_3Gene_Sim[c(1:21, 32, 39:41)] %>%
  unique()%>% #remove duplicate values
  dplyr::select(c(5, 10, 11, 22:25))%>% #extract only columns needed for downstream analysis
  dplyr::rename(read_supp = Reads)

dGen <- Annot_Genion_Sim[c(1:27, 37, 44:46)]%>%
  unique()%>% #remove duplicate values
  dplyr::select(c(1, 11, 19, 28:31))%>% #extract only columns needed for downstream analysis
  dplyr::rename(read_supp = V5,
         fusion.gene.id = V1)

dFS <- Annot_FusionSeeker_Sim[c(1:23, 33, 40:42)]%>%
  unique()%>% 
  dplyr::select(c(8, 15, 5, 24, 25:27))%>%
  dplyr::rename(read_supp = NumSupp)

dCTAT_LR<- Annot_CTATLR_Sim[c(1:25, 43, 50:52)] %>% 
  unique()%>% #remove duplicate values
  dplyr::select(c(2, 22, 23, 26:29))%>%#select only columns needed for downstream analysis
  dplyr::rename(read_supp = num_LR)

combined_data <- rbind(dLGF, dJF, dJF3, dGen, dFS, dCTAT_LR)
combined_data$read_supp <- as.numeric(combined_data$read_supp)

write_tsv(combined_data, paste0(Analysis_output_folder, "combined_data_FusionCalls.tsv"))

combined_data_annotation <- filter(combined_data, read_supp >= 2) %>% 
  unique() %>% mutate(recall_category = case_when(
  grepl("false", fusionType) ~ "False_Call",
  grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType) ~ "Partial_Recall",
  TRUE ~ "True_Recall"))%>%
  group_by(control, depth, Sequence_Identity, Algorithm, recall_category) %>%
  summarise(RECALL_COUNT = n()) %>%
  summarise(False_Call_Number = sum(RECALL_COUNT[recall_category == "False_Call"], na.rm = TRUE),
            Partial_Call_Number = sum(RECALL_COUNT[recall_category == "Partial_Recall"], na.rm = TRUE),
            True_Call_Number = sum(RECALL_COUNT[recall_category == "True_Recall"], na.rm = TRUE))%>% 
  ungroup() %>% 
  complete(control, depth, Sequence_Identity, Algorithm, fill = list(False_Call_Number = 0, Partial_Call_Number = 0, True_Call_Number = 0))
combined_data_annotation$depth <- factor(combined_data_annotation$depth, levels = c("1GB", "10GB", "100GB"))
combined_data_annotation$Sequence_Identity <- factor(combined_data_annotation$Sequence_Identity, levels = c("85%", "90%", "95%"))

write_tsv(combined_data_annotation, file=paste0(Analysis_output_folder, "combined_data_annotation.tsv"))
