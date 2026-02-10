##############################
# Fusion Types Overlap Upset #
##############################
FUSIONLIST <- c("JAFFAL", "Genion", "LongGF", "FusionSeeker", "CTAT-LR-Fusion", "GFSeeker")
## Make Supplementary Figure 3 
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  GFSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "GFSeeker", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("false", fusionType), read_supp >= 2)$fusion.gene.id))

Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)

# Convert list to binary membership matrix
Fusion_overlap_matrix <- as.data.frame(Fusion_overlap)  # Convert list-based to a dataframe

# Ensure the matrix is in the correct format (0/1 values)
Fusion_overlap_matrix <- as.data.frame(lapply(Fusion_overlap_matrix, as.integer))

# Generate the UpSet plot using ComplexUpset
False_Overlap <- ComplexUpset::upset(
  Fusion_overlap_matrix, names(Fusion_overlap_matrix)[1:num_ALG],  # First four sets
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=FUSIONLIST, color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='red'), guide='none') 
    + ylab('False Calls')
    )))

#Supplementary Figure 3B
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  GFSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "GFSeeker", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", grepl("truncated|reverse|chromosomal_misalignment|wrong", fusionType), read_supp >= 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)
# Generate an upset plot without specifying intersect, which uses the default behavior
PartialCall_Overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:num_ALG]),
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=FUSIONLIST, color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='black'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='orange'), guide='none')
    + ylab('Partial Recall')
    ))) 
supp_Figure3 <-(False_Overlap/PartialCall_Overlap)+   
  plot_annotation(tag_levels = list(c("A", "", "", "B")))
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure3.pdf", plot = supp_Figure3, width = 297, height =210 , units = "mm")
ggsave("~/LongReadFusionCallerBenchmark/Figures/supp_Figure3.png", plot = supp_Figure3, width = 297, height =210 , units = "mm")