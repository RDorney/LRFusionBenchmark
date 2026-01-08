######################
# True Calls Overlap #
######################
True_overlap <- ComplexUpset::upset(
  Fusion_overlap, intersect = names(Fusion_Simulated_Fusion_overlap_sets[1:4]), name='',
  width_ratio=0.1, sort_intersections_by='degree',
  queries=list(upset_query(intersect=c("JAFFAL", "Genion", "LongGF", "FusionSeeker"), 
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