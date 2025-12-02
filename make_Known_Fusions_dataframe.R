## Huh7 to Liver fusion Recall
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggbreak)
library(ComplexUpset)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
####Load biomaRT to gene Ensembl gene IDS####
library(biomaRt)
ensemblv113 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",  host = "https://oct2024.archive.ensembl.org")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

####Format Known Fusions####
#Collate known fusions
KF_ChimerKB <- ChimerDB_KB%>% filter(!is.na(Disease), Disease!="NA") %>% 
  select(c(Fusion_pair, H_gene, T_gene, Source, Disease, Validation, webSource, database, subdatabase, order)) %>% unique() %>%
  rename(known_fusion = Fusion_pair, Gene1 = H_gene, Gene2 = T_gene)%>%
  mutate(known_fusion = gsub("-", "::", known_fusion))
KF_ChimerSeq4 <- ChimerDB_Seq4 %>% filter(!is.na(Cancertype), Cancertype!="NA")  %>%
  select(c(Fusion_pair, H_gene, T_gene, Cancertype, Source, webSource, database, subdatabase, order)) %>%
  unique() %>%
  rename(known_fusion = Fusion_pair, Gene1 = H_gene, Gene2 = T_gene, Disease = Cancertype) %>%
  mutate(known_fusion = gsub("-", "::", known_fusion)) 
KF_ChimerPub4 <- ChimerDB_Pub4 %>% filter(!is.na(Disease), Disease!="NA") %>% 
  select(c(Fusion_pair, H_gene, T_gene, Disease, Validation, database, subdatabase, order)) %>%
  unique() %>%
  rename(known_fusion = Fusion_pair, Gene1 = H_gene, Gene2 = T_gene) %>%
  mutate(known_fusion = gsub("-", "::", known_fusion)) 

KF_Cosmic <- Cosmic_Fusion %>% 
  select(c(FIVE_PRIME_GENE_SYMBOL, FIVE_PRIME_TRANSCRIPT_ID, THREE_PRIME_GENE_SYMBOL, THREE_PRIME_TRANSCRIPT_ID, database, subdatabase, order)) %>% unique() %>%
  rename(Gene1 = FIVE_PRIME_GENE_SYMBOL, Gene2 = THREE_PRIME_GENE_SYMBOL) %>%
  unite("known_fusion" , Gene1:Gene2, sep = "::", remove = FALSE)

KF_DepMap_Liv <- DepMap_Liver %>%
  extract(Gene1, into = c("Gene1_name", "Gene1_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  extract(Gene2, into = c("Gene2_name", "Gene2_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  unite("fusion_name", Gene1_name,Gene2_name, sep="::", remove = FALSE ) %>% 
  select(c(fusion_name, Gene1_name, Gene1_ensembl, Gene2_name, Gene2_ensembl, database, subdatabase, order))%>%
  unique()%>%
  rename(known_fusion = fusion_name)
KF_DepMap_HCC <- DepMap_HCC %>%
  extract(Gene1, into = c("Gene1_name", "Gene1_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  extract(Gene2, into = c("Gene2_name", "Gene2_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  unite("fusion_name", Gene1_name,Gene2_name, sep="::", remove = FALSE ) %>% 
  select(c(fusion_name, Gene1_name, Gene1_ensembl, Gene2_name, Gene2_ensembl, database, subdatabase, order)) %>%
  unique()%>%
  rename(known_fusion = fusion_name, Gene1 = Gene1_name, Gene2 = Gene2_name)
KF_DepMap_Huh7 <- DepMap_Huh7 %>%
  extract(Gene1, into = c("Gene1_name", "Gene1_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  extract(Gene2, into = c("Gene2_name", "Gene2_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  unite("fusion_name", Gene1_name,Gene2_name, sep="::", remove = FALSE ) %>% 
  select(c(fusion_name, Gene1_name, Gene1_ensembl, Gene2_name, Gene2_ensembl, database, subdatabase, order)) %>%
  unique()%>%
  rename(known_fusion = fusion_name, Gene1 = Gene1_name, Gene2 = Gene2_name)

KF_FPIA <- FPIA_tcga_fusion %>% filter(!is.na(cancer) ,  cancer!="NA") %>% 
  separate(fusion, into = c("fusion_genea", "fusion_geneb"), sep = "--", remove = FALSE) %>%
  select(c(fusion, cancer, fusion_genea, fusion_geneb, database, subdatabase, order))%>%
  unique()%>% rename(known_fusion = fusion, Gene1 = fusion_genea, Gene2 = fusion_geneb, Disease = cancer) %>%
  mutate(known_fusion = gsub("--", "::", known_fusion)) 
KF_Mitelman <- mitelmandb_genefusion %>% 
  select(c("Fusion Gene", "Tumor Morphology", database, subdatabase, order))%>%
  separate(col = "Fusion Gene", into = c("Gene1", "Gene2"), sep ="::", remove = FALSE ) %>%
  unique()%>%
  rename(known_fusion = "Fusion Gene", Disease = "Tumor Morphology")

Known_Fusions_df <- full_join(KF_ChimerKB, KF_ChimerSeq4) %>%
  full_join(KF_ChimerPub4) %>% full_join(KF_Cosmic) %>%
  full_join(KF_DepMap_Liv) %>% full_join(KF_DepMap_HCC) %>%
  full_join(KF_DepMap_Huh7) %>% full_join(KF_FPIA) %>%
  full_join(KF_Mitelman)

###Search for gene IDs based on gene name
gene_info_KnownFusions <- rbind(getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', "external_synonym", "chromosome_name"),
  filters = 'external_gene_name',
  values = unique(c(Known_Fusions_df$Gene1, Known_Fusions_df$Gene2)),
  mart = ensembl
), getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', "external_synonym", "chromosome_name"),
  filters = "external_synonym",
  values = unique(c(Known_Fusions_df$Gene1, Known_Fusions_df$Gene2)),
  mart = ensembl
),getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', "external_synonym", "chromosome_name"),
  filters = 'external_gene_name',
  values = unique(c(Known_Fusions_df$Gene1, Known_Fusions_df$Gene2)),
  mart = ensemblv113
), getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', "external_synonym", "chromosome_name"),
  filters = "external_synonym",
  values = unique(c(Known_Fusions_df$Gene1, Known_Fusions_df$Gene2)),
  mart = ensemblv113
))%>% filter(!str_detect(chromosome_name, "PATCH|HSCHR")) %>% pivot_longer(cols = external_gene_name:external_synonym, names_to = "name_type", values_to = "gene_names")%>%
  unique()



