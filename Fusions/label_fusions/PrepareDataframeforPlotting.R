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
  mutate(spanning.reads = split_reads1 + split_reads2)

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

