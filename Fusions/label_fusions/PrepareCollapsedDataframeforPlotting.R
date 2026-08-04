#############################################
# Prepare Dataframe for plotting (Collapsed)
#############################################
Annot_CTATLR_Huh7_Discovery_collapseddf <- Annot_CTATLR_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_LongGF_Huh7_Discovery_collapseddf <- Annot_LongGF_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_Genion_Huh7_Discovery_collapseddf<- Annot_Genion_Huh7_Discovery_collapsed
Annot_Genion_Huh7_Discovery_collapseddf$fusion.ens.gene.id <- Annot_Genion_Huh7_Discovery_collapseddf$V1

Annot_FusionSeeker_Huh7_Discovery_collapseddf <- Annot_FusionSeeker_Huh7_Discovery_collapsed
Annot_FusionSeeker_Huh7_Discovery_collapsed$fusion.ens.gene.id <- Annot_FusionSeeker_Huh7_Discovery_collapseddf$fusion.gene.id

Annot_GFSeeker_Huh7_Discovery_collapseddf <- Annot_GFSeeker_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_JAFFAL_Huh7_Discovery_collapseddf <- Annot_JAFFAL_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

JAFFAL3Gene_Huh7_Discovery_collapseddf <- JAFFAL3Gene_Huh7_Discovery_collapsed %>%
  unite("fusion.ens.gene.id", ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id, 
        sep = "::", remove = FALSE)

Annot_STARFusion_Huh7_Discovery_collapseddf <- Annot_STARFusion_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(GENEID1, "::", GENEID2))

Annot_Arriba_Huh7_Discovery_collapseddf <- Annot_Arriba_Huh7_Discovery_collapsed %>%
  mutate(fusion.ens.gene.id = paste0(GENEID1, "::", GENEID2))

cols <- c("library_type", "Platform", 
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id")
df_list <- list(
  Annot_CTATLR_Huh7_Discovery_collapseddf, Annot_Genion_Huh7_Discovery_collapseddf
  , Annot_LongGF_Huh7_Discovery_collapseddf, 
  Annot_FusionSeeker_Huh7_Discovery_collapseddf, Annot_GFSeeker_Huh7_Discovery_collapseddf, 
  Annot_JAFFAL_Huh7_Discovery_collapseddf, JAFFAL3Gene_Huh7_Discovery_collapseddf,
  Annot_STARFusion_Huh7_Discovery_collapseddf, Annot_Arriba_Huh7_Discovery_collapseddf
)

fusion_Huh7_discovery_collapsed <- df_list %>%
  map(~ .x %>% unique() %>% select(all_of(cols))) %>%
  bind_rows()

write_tsv(fusions_Huh7_discovery_collapsed, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/fusions_Huh7_discovery.tsv")
write_tsv(fusions_Huh7_discovery_collapsed, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_Huh7_discovery.tsv.gz")

CTATLR_Huh7_Discovery_collapsed <- Annot_CTATLR_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = num_LR)
CTATLR_Huh7_Discovery_collapsed$spanning.pairs <- NA

LongGF_Huh7_Discovery_collapsed <- Annot_LongGF_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = V3)
LongGF_Huh7_Discovery_collapsed$spanning.pairs <- NA

Genion_Huh7_Discovery_collapsed <- Annot_Genion_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = V5)
Genion_Huh7_Discovery_collapsed$spanning.pairs <- NA

FusionSeeker_Huh7_Discovery_collapsed <- Annot_FusionSeeker_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = NumSupp)
FusionSeeker_Huh7_Discovery_collapsed$spanning.pairs <- NA

GFSeeker_Huh7_Discovery_collapsed <- Annot_GFSeeker_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = "support num")
GFSeeker_Huh7_Discovery_collapsed$spanning.pairs <- NA

JAFFAL_Huh7_Discovery_collapsed <- Annot_JAFFAL_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = spanning.reads) 

JAFFAL3Gene_Huh7_Discovery_collapsed <- JAFFAL3Gene_Huh7_Discovery_collapseddf %>%
  rename(spanning.reads = Reads) 
JAFFAL3Gene_Huh7_Discovery$spanning.pairs <- NA

cols <- c("library_type", "Platform",
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id", 
          "spanning.reads", "spanning.pairs")

df_list <- list(
  CTATLR_Huh7_Discovery_collapsed, Genion_Huh7_Discovery_collapsed, 
  LongGF_Huh7_Discovery_collapsed, FusionSeeker_Huh7_Discovery_collapsed, 
  GFSeeker_Huh7_Discovery_collapsed, JAFFAL_Huh7_Discovery_collapsed, JAFFAL3Gene_Huh7_Discovery_collapsed
)

fusions_rs_Huh7_discovery_collapsed <- df_list %>%
  map(~ .x %>% unique() %>% select(all_of(cols))) %>%
  bind_rows()

write_tsv(fusions_rs_Huh7_discovery_collapsed, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_readsupport_Huh7_discovery.tsv.gz")

#########################
# Plotting
#########################
ggplot(dplyr::filter(fusions_Huh7_discovery, !fusionType %in% c(
  "Mitochondrial:Genomic",
  "Mitochondrial:Mitochondrial",
  "Self-Misalignment"
)), 
aes(x = Discovery, colour = Algorithm, shape=RNA_sample)) +
  geom_point(stat = "count") +
  theme_bw() +
  labs(title = "", subtitle = "minimum read support of 2", x = "", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(library_type~ Platform, scales = "free_y") +
  labs(fill = "Algorithm")+
  scale_colour_manual(values = Alg_colour_map)+
  scale_y_continuous(labels = label_comma()) 

count_data_discovery <- fusions_Huh7_discovery %>%
  filter(!fusionType %in% c(
    "Mitochondrial:Genomic",
    "Mitochondrial:Mitochondrial",
    "Self-Misalignment"
  )) %>%
  count(RNA_sample, Algorithm, fusionType, library_type, Platform, Discovery )

ggplot(
  count_data_discovery,
  aes(x = RNA_sample, y = n)
) +
  geom_boxplot(aes(group = RNA_sample), outlier.shape = NA) +
  geom_jitter(aes(colour = Algorithm), width = 0.2, size = 2) +
  theme_bw() +
  labs(
    title = "",
    subtitle = "minimum read support of 2",
    x = "",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  facet_grid(fusionType ~ library_type + Platform) +
  scale_y_log10() +
  scale_colour_manual(values = Alg_colour_map)


count_data_discovery <- fusions_Huh7_discovery %>%
  filter(!fusionType %in% c(
    "Mitochondrial:Genomic",
    "Mitochondrial:Mitochondrial",
    "Self-Misalignment"
  )) %>%
  count(RNA_sample, Algorithm, library_type, Platform, Discovery, Sample, Validation )


ggplot(
  count_data_discovery,
  aes(x = Platform, y = n)
) +
  geom_boxplot(aes(group = Platform), outlier.shape = NA) +
  geom_jitter(aes(colour = Algorithm), width = 0.2, size = 2) +
  theme_bw() +
  labs(
    title = "",
    subtitle = "minimum read support of 2",
    x = "",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  facet_grid(Sample+Discovery~library_type, scales = "free_y") +
  #scale_y_log10() +
  scale_colour_manual(values = Alg_colour_map)

ggplot(
  count_data_discovery,
  aes(x = Platform, y = n)
) +
  geom_boxplot(aes(group = Platform), outlier.shape = NA) +
  geom_jitter(aes(colour = Algorithm), width = 0.2, size = 2) +
  theme_bw() +
  labs(
    title = "",
    subtitle = "minimum read support of 2",
    x = "",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  facet_grid(Sample+Discovery~library_type, scales = "free_y") +
  scale_y_log10() +
  scale_colour_manual(values = Alg_colour_map)

