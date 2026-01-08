#########################
# Generate overlap sets #
#########################
# Generate upset plots to see the total number of overlapping fusions between fusion callers
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)