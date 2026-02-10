#Append Depth & SeqID data

Annot_LongGF_Sim <- Annot_LongGF_Sim %>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB" # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_JAFFAL_Sim <- Annot_JAFFAL_Sim%>%
  mutate(depth = case_when(
    grepl("100GB", sample) ~ "100GB",
    grepl("10GB", sample) ~ "10GB",
    grepl("1GB", sample) ~ "1GB"# Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", sample) ~ "95%",
    grepl("90", sample) ~ "90%",
    grepl("85", sample) ~ "85%" # Keep original value if no match
  ))
Annot_JAFFAL_3Gene_Sim <- Annot_JAFFAL_3Gene_Sim%>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"# Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_Genion_Sim <- Annot_Genion_Sim %>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_FusionSeeker_Sim <- Annot_FusionSeeker_Sim %>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_GFSeeker_Sim <- Annot_GFSeeker_Sim %>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_CTATLR_Sim <- Annot_CTATLR_Sim %>%
  mutate(depth = case_when(
    grepl("100GB", Source) ~ "100GB",
    grepl("10GB", Source) ~ "10GB",
    grepl("1GB", Source) ~ "1GB"  # Keep original value if no match
  ), Sequence_Identity = case_when(
    grepl("95", Source) ~ "95%",
    grepl("90", Source) ~ "90%",
    grepl("85", Source) ~ "85%" # Keep original value if no match
  ))
Annot_LongGF_Sim$Algorithm <- "LongGF"
Annot_JAFFAL_Sim$Algorithm <- "JAFFAL"
Annot_JAFFAL_3Gene_Sim$Algorithm <- "JAFFAL"
Annot_Genion_Sim$Algorithm <- "Genion"
Annot_FusionSeeker_Sim$Algorithm <- "FusionSeeker"
Annot_CTATLR_Sim$Algorithm <- "CTAT-LR-Fusion"
Annot_GFSeeker_Sim$Algorithm <- "GFSeeker"