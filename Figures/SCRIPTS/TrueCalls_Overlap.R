######################
# True Calls Overlap #
######################
# Generate upset plots to see the total number of overlapping fusions between fusion callers
Fusion_Simulated_Fusion_overlap_sets <- list(
  JAFFAL = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "JAFFAL", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  Genion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "Genion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  LongGF = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "LongGF", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  FusionSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "FusionSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  GFSeeker = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "GFSeeker", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id),
  CTAT_LR_Fusion = unique(filter(combined_data, depth == "100GB", Sequence_Identity == "95%", control == "positive", Algorithm == "CTAT-LR-Fusion", !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType), read_supp >= 2)$fusion.gene.id)
)

# Convert the list to a format compatible with UpSetR
Fusion_overlap <- fromList(Fusion_Simulated_Fusion_overlap_sets)

True_overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:num_ALG]), name='',
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker", "GFSeeker"), 
                           color='blue')),
  base_annotations = list(
    'Intersection size'=(intersection_size(
      text_colors=c(on_background='black', on_bar='white'),
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='springgreen4'), guide='none')
    + ylab('True Recall')
    ))) 
ggsave("~/LongReadFusionCallerBenchmark/Figures/Figure1H.pdf", 
       plot = True_overlap, width = 297, height = 210, units = "mm")
