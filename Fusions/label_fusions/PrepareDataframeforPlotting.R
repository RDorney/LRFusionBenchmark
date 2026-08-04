#################################
# Prepare Dataframe for plotting
#################################
library(purrr)
library(dplyr)
library(tidyr)
Annot_CTATLR_Huh7_Discoverydf <- Annot_CTATLR_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_LongGF_Huh7_Discoverydf <- Annot_LongGF_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_Genion_Huh7_Discoverydf <- Annot_Genion_Huh7_Discovery
Annot_Genion_Huh7_Discoverydf$fusion.ens.gene.id <- Annot_Genion_Huh7_Discoverydf$V1

Annot_FusionSeeker_Huh7_Discoverydf <- Annot_FusionSeeker_Huh7_Discovery
Annot_FusionSeeker_Huh7_Discoverydf$fusion.ens.gene.id <- Annot_FusionSeeker_Huh7_Discoverydf$fusion.gene.id

Annot_GFSeeker_Huh7_Discoverydf <- Annot_GFSeeker_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_JAFFAL_Huh7_Discoverydf <- Annot_JAFFAL_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

JAFFAL3Gene_Huh7_Discoverydf <- JAFFAL3Gene_Huh7_Discovery %>%
  unite("fusion.ens.gene.id", ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id, 
        sep = "::", remove = FALSE)

Annot_STARFusion_Huh7_Discoverydf <- Annot_STARFusion_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(GENEID1, "::", GENEID2))

Annot_Arriba_Huh7_Discoverydf <- Annot_Arriba_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(GENEID1, "::", GENEID2))

cols <- c("library_type", "Platform", 
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id")
df_list <- list(
  Annot_CTATLR_Huh7_Discoverydf, Annot_Genion_Huh7_Discoverydf, 
  Annot_LongGF_Huh7_Discoverydf, 
  Annot_FusionSeeker_Huh7_Discoverydf, Annot_GFSeeker_Huh7_Discoverydf, 
  Annot_JAFFAL_Huh7_Discoverydf, JAFFAL3Gene_Huh7_Discoverydf,
  Annot_STARFusion_Huh7_Discoverydf, Annot_Arriba_Huh7_Discoverydf
)

fusions_Huh7_discovery <- df_list %>%
  map(~ .x %>% unique() %>% 
        dplyr::select(all_of(cols))) %>%
  bind_rows()

write_tsv(fusions_Huh7_discovery, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/fusions_Huh7_discovery.tsv")
write_tsv(fusions_Huh7_discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_Huh7_discovery.tsv.gz")

CTATLR_Huh7_Discovery <- Annot_CTATLR_Huh7_Discoverydf%>%
  dplyr::rename(spanning.reads = num_LR)
CTATLR_Huh7_Discovery$spanning.pairs <- NA

LongGF_Huh7_Discovery <- Annot_LongGF_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = V3)
LongGF_Huh7_Discovery$spanning.pairs <- NA

Genion_Huh7_Discovery <- Annot_Genion_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = V5)
Genion_Huh7_Discovery$spanning.pairs <- NA

FusionSeeker_Huh7_Discovery <- Annot_FusionSeeker_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = NumSupp)
FusionSeeker_Huh7_Discovery$spanning.pairs <- NA

GFSeeker_Huh7_Discovery <- Annot_GFSeeker_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = "support num")
GFSeeker_Huh7_Discovery$spanning.pairs <- NA

JAFFAL_Huh7_Discovery <- Annot_JAFFAL_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = spanning.reads) 

JAFFAL3Gene_Huh7_Discovery <- JAFFAL3Gene_Huh7_Discoverydf %>%
  dplyr::rename(spanning.reads = Reads) 
JAFFAL3Gene_Huh7_Discovery$spanning.pairs <- NA

STARFusion_Huh7_Discovery <- Annot_STARFusion_Huh7_Discoverydf %>%
  dplyr::rename(spanning.pairs = SpanningFragCount) %>%
  dplyr::rename(spanning.reads = JunctionReadCount)

Arriba_Huh7_Discovery <- Annot_Arriba_Huh7_Discoverydf %>%
  dplyr::rename(spanning.pairs = discordant_mates) %>%
  mutate(spanning.reads = sum((split_reads1 + split_reads2)))

cols <- c("library_type", "Platform", 
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id", 
          "spanning.reads", "spanning.pairs")

df_list <- list(
  Arriba_Huh7_Discovery, STARFusion_Huh7_Discovery, 
  CTATLR_Huh7_Discovery, Genion_Huh7_Discovery, 
  LongGF_Huh7_Discovery, FusionSeeker_Huh7_Discovery, 
  GFSeeker_Huh7_Discovery, JAFFAL_Huh7_Discovery, JAFFAL3Gene_Huh7_Discovery
)

fusions_rs_Huh7_discovery <- df_list %>%
  map(~ .x %>% unique() %>% dplyr::select(all_of(cols))) %>%
  bind_rows()

write_tsv(fusions_rs_Huh7_discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_readsupport_Huh7_discovery.tsv.gz")

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

