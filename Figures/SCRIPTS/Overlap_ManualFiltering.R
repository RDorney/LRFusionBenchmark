#######################################
# Manual Filtering & Overlap Analysis #
#######################################
#Repeat analysis with removed obvious false chimeras to see if they have any influence. 
#Any differences observed are negligible
combined_data_with_heuristic <- filter(combined_data, read_supp >= 2) %>%  
  mutate(stat_category = case_when(
    grepl("false_fusion:mitochondrial|false_fusion:self_misalignment", fusionType) ~ "obvious_false",
    TRUE ~ "other")) %>% filter(stat_category == "other") %>% select(-'stat_category')

Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
# Generate an upset plot without specifying intersect, which uses the default behavior
True_overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:4]), name='',
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='white'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='springgreen4'), guide='none')
    + ylab('True Recall')
    ))) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/Compare_Contrast_True_Calls.pdf", plot = True_overlap, width = 297, height = 210, units = "mm")

Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
# Generate an upset plot without specifying intersect, which uses the default behavior
PartialCall_Overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:4]),
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='orange'), guide='none')
    + ylab('Partial Recall')
    ))) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/Compare_Contrast_Partial_Calls.pdf", plot = PartialCall_Overlap, width = 297, height = 210, units = "mm")

Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data_with_heuristic, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id)
)

Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
# Convert list to binary membership matrix
Fusion_overlap_matrix <- as.data.frame(Fusion_overlap)  # Convert list-based to a dataframe

# Ensure the matrix is in the correct format (0/1 values)
Fusion_overlap_matrix <- as.data.frame(lapply(Fusion_overlap_matrix, as.integer))

# Generate the UpSet plot using ComplexUpset
False_Overlap <- ComplexUpset::upset(
  Fusion_overlap_matrix, names(Fusion_overlap_matrix)[1:4],  # First four sets
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='red'), guide='none') 
    + ylab('False Calls')
    )))
ggsave("~/LongReadFusionCallerBenchmark/Figures/Compare_Contrast_False_Calls.pdf", plot = False_Overlap, width = 297, height = 210, units = "mm")
