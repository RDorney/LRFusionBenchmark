################################
# Title: Check Known Fusions
# Author: Ryley Dorney
# Date: Mar 2026
################################
library(tidyverse)
# Load in Known Fusions List
#source("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/ImportFormat_Known_Huh7.R")
huh7fusions_manualannot<- read.csv("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/KnownHuh7Fusion_ah.csv")

View(huh7fusions_manualannot)
huh7fusions_manualannot<- huh7fusions_manualannot%>% mutate(fusion.gene.id = paste0(GENEID.x, "::", GENEID.y))
#huh7fusions_manualannot$Discovery <- "Known"
#######################
# Check CTAT-LR-Fusion
#######################
Annot_CTATLR_Huh7_Known<- inner_join(Annot_CTATLR_Huh7, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "Known")

Annot_CTATLR_Huh7_Known<- right_join(Annot_CTATLR_Huh7_Known, Annot_CTATLR_Huh7)
idx <- which(is.na(Annot_CTATLR_Huh7_Known$Discovery) | Annot_CTATLR_Huh7_Known$Discovery == "")# Check if the discovery is empty

  # Reverse match
Annot_CTATLR_Huh7_RevKnown <-  Annot_CTATLR_Huh7_Known[idx,] %>%
  select(where(~ !all(is.na(.x)))) %>% 
  inner_join(huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
    mutate(Discovery = "Reverse Known") 

Annot_CTATLR_Huh7_Discovery <- left_join(Annot_CTATLR_Huh7, Annot_CTATLR_Huh7_Known) %>% 
  full_join(Annot_CTATLR_Huh7_RevKnown) %>%
  unique()
nrow(Annot_CTATLR_Huh7_Discovery)
idx <- which(is.na(Annot_CTATLR_Huh7_Discovery$Discovery) | Annot_CTATLR_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_CTATLR_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
unique(Annot_CTATLR_Huh7) %>% nrow()

write_tsv(Annot_CTATLR_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_CTATLR_Huh7_Discovery.tsv.gz")

Annot_CTATLR_Huh7_Discovery_collapsed <- Annot_CTATLR_Huh7_Discovery%>%
  dplyr::group_by(across(c("#FusionName","LeftGene","RightGene",
                           "chrom1","chrom2",
                           "Source", "Cell_Line","Algorithm",
                           "RNA_sample","library_type","Platform","fusionType",
                           "full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_read = sum(num_LR), .groups = "drop") 
write_tsv(Annot_CTATLR_Huh7_Discovery_collapsed, 
          file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_CTATLR_Huh7_Discovery_collapsed.tsv.gz")

#######################
# Check JAFFA/L
#######################
Annot_JAFFAL_Huh7_Known <-inner_join(unique(Annot_JAFFAL_Huh7), huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "Known")

Annot_JAFFAL_Huh7_Known<- right_join(Annot_JAFFAL_Huh7_Known, unique(Annot_JAFFAL_Huh7), 
                                     by=c("Gene1.x" = "Gene1", "Gene2.x" = "Gene2",
                                          "sample","fusion.genes","chrom1","base1","strand1","chrom2","base2","strand2",
                                          "gap..kb.","spanning.pairs","spanning.reads","inframe","aligns",
                                          "rearrangement","contig","contig.break","classification","known",
                                          "Algorithm","Platform","library_type","Cell_Line","RNA_sample",
                                          "full_label","ensembl_gene_id.x","G1START","G1END","ensembl_gene_id.y",
                                          "G2START","G2END","fusionType"))%>%
  unique()%>% rename("Gene1"="Gene1.x","Gene2"="Gene2.x")
nrow(unique(Annot_JAFFAL_Huh7))
nrow(Annot_JAFFAL_Huh7_Known)
idx <- which(is.na(Annot_JAFFAL_Huh7_Known$Discovery) | Annot_JAFFAL_Huh7_Known$Discovery == "")# Check if the discovery is empty

# Reverse match
Annot_JAFFAL_Huh7_RevKnown <-  Annot_JAFFAL_Huh7_Known[idx,] %>%
  select(where(~ !all(is.na(.x)))) %>%
  inner_join(huh7fusions_manualannot, 
              by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
  rename("Gene1"="Gene1.x","Gene2"="Gene2.x")%>% 
  mutate(Discovery = "Reverse Known") 
idx <- !which(is.na(Annot_JAFFAL_Huh7_Known$Discovery) | Annot_JAFFAL_Huh7_Known$Discovery == "")# Check if the discovery is empty

Annot_JAFFAL_Huh7_Discovery <- rbind(Annot_JAFFAL_Huh7_RevKnown, Annot_JAFFAL_Huh7_Known) %>%
  filter(!is.na(Discovery)) %>% right_join(Annot_JAFFAL_Huh7)%>%
  unique()
nrow(Annot_JAFFAL_Huh7_Discovery)

idx <- which(is.na(Annot_JAFFAL_Huh7_Discovery$Discovery) | Annot_JAFFAL_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_JAFFAL_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
Annot_JAFFAL_Huh7_Discovery %>% nrow()

#write_tsv(Annot_JAFFAL_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_JAFFAL_Huh7_Discovery.tsv.gz")
write_tsv(Annot_JAFFAL_Huh7_Discovery, xzfile("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_JAFFAL_Huh7_Discovery.tsv.xz", compression = 9))

Annot_JAFFAL_Huh7_Discovery_collasped <- Annot_JAFFAL_Huh7_Discovery %>%
  dplyr::group_by(across(c("sample","fusion.genes","Gene1","Gene2",
                           "chrom1","chrom2",
                           "Cell_Line","Algorithm", "fusionType",
                           "RNA_sample","library_type","Platform",
                           "full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_pairs = sum(spanning.pairs), tot_span_read = sum(spanning.reads), .groups = "drop") 
write_tsv(Annot_JAFFAL_Huh7_Discovery_collasped, 
          xzfile("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_JAFFAL_Huh7_Discovery_collasped.tsv.xz", compression = 9))

#check 3 gene
JAFFAL3Gene_Huh7_12 <- inner_join(Annot_JAFFAL_Huh7_3Gene, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "contains Known 1.2")

JAFFAL3Gene_Huh7_23 <- Annot_JAFFAL_Huh7_3Gene %>%
  inner_join(huh7fusions_manualannot, 
            by = c("ensembl_gene_id.y" = "GENEID.x", "ensembl_gene_id" = "GENEID.y")) %>%
    mutate(Discovery = "contains Known 2.3")%>%
  select(where(~ !all(is.na(.x))))

JAFFAL3Gene_Huh7_Discovery<-rbind(JAFFAL3Gene_Huh7_12, JAFFAL3Gene_Huh7_23)%>%
  unique()%>%
  rename("Gene1"="Gene1.x","Gene2"="Gene2.x")%>% 
  right_join(Annot_JAFFAL_Huh7_3Gene)%>%
  # If Discovery is still NA, it didn't find a match in either direction
  mutate(Discovery = if_else(is.na(Discovery), "Putative Novel", Discovery))

write_tsv(JAFFAL3Gene_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/JAFFAL3Gene_Huh7_Discovery.tsv.gz")

#######################
# Check LongGF
#######################
LongGF_Huh7_Known <- inner_join(unique(Annot_LongGF_Huh7), huh7fusions_manualannot, 
           by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
  mutate(Discovery = "Known")

LongGF_Huh7_Known_check<- right_join(LongGF_Huh7_Known, Annot_LongGF_Huh7)
idx <- which(is.na(LongGF_Huh7_Known_check$Discovery) | LongGF_Huh7_Known_check$Discovery == "")# Check if the discovery is empty

# Reverse match
LongGF_Huh7_RevKnown <-  LongGF_Huh7_Known_check[idx,] %>%
  select(where(~ !all(is.na(.x)))) %>% 
  inner_join(huh7fusions_manualannot, 
             by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
  mutate(Discovery = "Reverse Known") 

Annot_LongGF_Huh7_Discovery <- rbind(LongGF_Huh7_RevKnown, LongGF_Huh7_Known) %>%
  filter(!is.na(Discovery)) %>% right_join(Annot_LongGF_Huh7)%>%
  unique()
nrow(Annot_LongGF_Huh7_Discovery)

idx <- which(is.na(Annot_LongGF_Huh7_Discovery$Discovery) | Annot_LongGF_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_LongGF_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
Annot_LongGF_Huh7_Discovery %>% nrow()
nrow(Annot_LongGF_Huh7)

write_tsv(Annot_LongGF_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_LongGF_Huh7_Discovery.tsv.gz")

Annot_LongGF_Huh7_Discovery_collapsed <- Annot_LongGF_Huh7_Discovery %>%
  dplyr::group_by(across(c("V1","Source","V2",
                           "chromosome1","chromosome2",
                           "Cell_Line","Algorithm",
                           "RNA_sample","library_type","Platform",
                           "fusionType","full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_read = sum(V3), .groups = "drop")
write_tsv(Annot_LongGF_Huh7_Discovery_collapsed, 
          file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_LongGF_Huh7_Discovery_collapsed.tsv.gz")

#######################
# Check FusionSeeker
#######################
FusionSeeker_Huh7_Known <- inner_join(unique(Annot_FusionSeeker_Huh7), huh7fusions_manualannot, 
                                      by = c("Gene1" = "GENEID.x", "Gene2" = "GENEID.y")) %>%
  mutate(Discovery = "Known")
FusionSeeker_Huh7_Known_check<- right_join(FusionSeeker_Huh7_Known, Annot_FusionSeeker_Huh7)
idx <- which(is.na(FusionSeeker_Huh7_Known_check$Discovery) | FusionSeeker_Huh7_Known_check$Discovery == "")# Check if the discovery is empty

# Reverse match
FusionSeeker_Huh7_RevKnown <-  FusionSeeker_Huh7_Known_check[idx,] %>%
  select(where(~ !all(is.na(.x)))) %>% 
  inner_join(huh7fusions_manualannot, 
             by = c("Gene1" = "GENEID.y", "Gene2" = "GENEID.x")) %>%
  mutate(Discovery = "Reverse Known") 

Annot_FusionSeeker_Huh7_Discovery <- rbind(FusionSeeker_Huh7_RevKnown, FusionSeeker_Huh7_Known) %>%
  filter(!is.na(Discovery)) %>% right_join(Annot_FusionSeeker_Huh7)%>%
  unique()
nrow(Annot_FusionSeeker_Huh7_Discovery)

idx <- which(is.na(Annot_FusionSeeker_Huh7_Discovery$Discovery) | Annot_FusionSeeker_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_FusionSeeker_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
Annot_FusionSeeker_Huh7_Discovery %>% nrow()
nrow(Annot_FusionSeeker_Huh7)

write_tsv(Annot_FusionSeeker_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_FusionSeeker_Huh7_Discovery.tsv.gz")

Annot_FusionSeeker_Huh7_Discovery_collapsed <- Annot_FusionSeeker_Huh7_Discovery %>%
  dplyr::group_by(across(c("Source","fusionGene",
                           "Chrom1","Chrom2",
                           "Cell_Line","Algorithm",
                           "RNA_sample","library_type","Platform",
                           "fusionType","full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_read = sum(NumSupp), .groups = "drop") 

write_tsv(Annot_FusionSeeker_Huh7_Discovery_collapsed, 
          file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_FusionSeeker_Huh7_Discovery_collapsed.tsv.gz")


#######################
# Check GFSeeker
#######################
GFSeeker_Huh7_Known <- inner_join(unique(Annot_GFSeeker_Huh7), huh7fusions_manualannot, 
              by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
  mutate(Discovery = "Known")

GFSeeker_Huh7_Known_check<- right_join(GFSeeker_Huh7_Known, Annot_GFSeeker_Huh7)
idx <- which(is.na(GFSeeker_Huh7_Known_check$Discovery) | GFSeeker_Huh7_Known_check$Discovery == "")# Check if the discovery is empty

# Reverse match
GFSeeker_Huh7_RevKnown <-  GFSeeker_Huh7_Known_check[idx,] %>%
  select(where(~ !all(is.na(.x)))) %>% 
  inner_join(huh7fusions_manualannot, 
             by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
  mutate(Discovery = "Reverse Known") 

Annot_GFSeeker_Huh7_Discovery <- rbind(GFSeeker_Huh7_RevKnown, GFSeeker_Huh7_Known) %>%
  filter(!is.na(Discovery)) %>% right_join(Annot_GFSeeker_Huh7)%>%
  unique()
nrow(Annot_GFSeeker_Huh7_Discovery)

idx <- which(is.na(Annot_GFSeeker_Huh7_Discovery$Discovery) | Annot_GFSeeker_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_GFSeeker_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
Annot_GFSeeker_Huh7_Discovery %>% nrow()
nrow(Annot_GFSeeker_Huh7)

write_tsv(Annot_GFSeeker_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_GFSeeker_Huh7_Discovery.tsv.gz")

Annot_GFSeeker_Huh7_Discovery_collapsed <- Annot_GFSeeker_Huh7_Discovery%>%
  dplyr::group_by(across(c("gene1_name","gene2_name",
                           "chrom1","chrom2",
                           "Source", "Cell_Line","Algorithm",
                           "RNA_sample","library_type","Platform", 
                           "fusionType","full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_read = sum(`support num`), .groups = "drop") 
write_tsv(Annot_GFSeeker_Huh7_Discovery_collapsed, 
          file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_GFSeeker_Huh7_Discovery_collapsed.tsv.gz")

#######################
# Check Genion
#######################
Genion_Huh7_Known <- inner_join(unique(filter(Annot_Genion_Huh7, chr3 == "")), huh7fusions_manualannot, 
                                      by = c("V1.1" = "GENEID.x", "V1.2" = "GENEID.y")) %>%
  mutate(Discovery = "Known")

Genion_Huh7_Known_check<- right_join(Genion_Huh7_Known, Annot_Genion_Huh7)
idx <- which(is.na(Genion_Huh7_Known_check$Discovery) | Genion_Huh7_Known_check$Discovery == "")# Check if the discovery is empty

# Reverse match
Genion_Huh7_RevKnown <-  Genion_Huh7_Known_check[idx,] %>% 
  unique() %>%
  select(where(~ !all(is.na(.x)))) %>% 
  filter(chr3 == "")%>%
  inner_join(huh7fusions_manualannot, 
             by = c("V1.1" = "GENEID.y", "V1.2" = "GENEID.x")) %>%
  mutate(Discovery = "Reverse Known") 

Genion_Huh7_check<- right_join(rbind(Genion_Huh7_Known, Genion_Huh7_RevKnown), Annot_Genion_Huh7)
idx <- which(is.na(Genion_Huh7_check$Discovery) | Genion_Huh7_check$Discovery == "")# Check if the discovery is empty

Genion_Huh7_multi_Known12  <-  Genion_Huh7_check[idx,] %>% 
  unique() %>%
  filter(chr3 != "")%>%
  select(where(~ !all(is.na(.x)))) %>%
  inner_join(huh7fusions_manualannot, 
             by = c("V1.1" = "GENEID.y", "V1.2" = "GENEID.x")) %>%
  mutate(Discovery = "contains Known 1.2") 

Genion_Huh7_multi_Known23 <-  Genion_Huh7_check[idx,] %>% 
  unique() %>%
  filter(chr3 != "")%>%
  select(where(~ !all(is.na(.x)))) %>% 
  inner_join(huh7fusions_manualannot, 
             by = c("V1.2" = "GENEID.y", "V1.3" = "GENEID.x")) %>%
  mutate(Discovery = "contains Known 2.3") 
Genion_Huh7_multi_Known34 <-  Genion_Huh7_check[idx,] %>% 
  unique() %>% 
  filter(chr4 != "")%>%
  select(where(~ !all(is.na(.x)))) %>%
  inner_join(huh7fusions_manualannot, 
             by = c("V1.3" = "GENEID.y", "V1.4" = "GENEID.x")) %>%
  mutate(Discovery = "contains Known 3.4")

# no recall in the multi gene fusions
Annot_Genion_Huh7_Discovery <- rbind(Genion_Huh7_RevKnown, Genion_Huh7_Known) %>%
  filter(!is.na(Discovery)) %>% right_join(Annot_Genion_Huh7)%>%
  unique()
nrow(Annot_Genion_Huh7_Discovery)

idx <- which(is.na(Annot_Genion_Huh7_Discovery$Discovery) | Annot_Genion_Huh7_Discovery$Discovery == "")# Check if the discovery is empty
Annot_Genion_Huh7_Discovery$Discovery[idx] <- "Putative Novel"
Annot_Genion_Huh7_Discovery %>% nrow()
nrow(Annot_Genion_Huh7)

write_tsv(Annot_Genion_Huh7_Discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_Genion_Huh7_Discovery.tsv.gz")

Annot_Genion_Huh7_Discovery_collapsed <- Annot_Genion_Huh7_Discovery%>%
  dplyr::group_by(across(c("Source","V1",
                           "V8",
                           "Cell_Line","Algorithm",
                           "RNA_sample","library_type","Platform",
                           "fusionType","full_label", "Discovery"
  ))) %>%
  dplyr::summarise(tot_span_read = sum(V5), .groups = "drop") 

write_tsv(Annot_Genion_Huh7_Discovery_collapsed, 
          file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Discovery/Annot_Genion_Huh7_Discovery_collapsed.tsv.gz")


#################################
# Prepare Dataframe for plotting
#################################
Annot_CTATLR_Huh7_Discovery <- Annot_CTATLR_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_LongGF_Huh7_Discovery <- Annot_LongGF_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_Genion_Huh7_Discovery$fusion.ens.gene.id <- Annot_Genion_Huh7_Discovery$V1
Annot_FusionSeeker_Huh7_Discovery$fusion.ens.gene.id <- Annot_FusionSeeker_Huh7_Discovery$fusion.gene.id.x

Annot_GFSeeker_Huh7_Discovery <- Annot_GFSeeker_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

Annot_JAFFAL_Huh7_Discovery <- Annot_JAFFAL_Huh7_Discovery %>%
  mutate(fusion.ens.gene.id = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y))

JAFFAL3Gene_Huh7_Discovery <- JAFFAL3Gene_Huh7_Discovery %>%
  unite("fusion.ens.gene.id", ensembl_gene_id.x, ensembl_gene_id.y, ensembl_gene_id, 
        sep = "::", remove = FALSE)

cols <- c("library_type", "Platform", "full_label", 
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id")
df_list <- list(
  Annot_CTATLR_Huh7_Discovery, Annot_Genion_Huh7_Discovery
  , Annot_LongGF_Huh7_Discovery, 
  Annot_FusionSeeker_Huh7_Discovery, Annot_GFSeeker_Huh7_Discovery, 
  Annot_JAFFAL_Huh7_Discovery, JAFFAL3Gene_Huh7_Discovery
)

fusion_Huh7_discovery <- df_list %>%
  map(~ .x %>% unique() %>% select(all_of(cols))) %>%
  bind_rows()

write_tsv(fusions_Huh7_discovery, file = "/bioinformatics/ryley/Gencode44/Huh7_Library/fusions_Huh7_discovery.tsv")
write_tsv(fusions_Huh7_discovery, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_Huh7_discovery.tsv.gz")

CTATLR_Huh7_Discovery <- Annot_CTATLR_Huh7_Discovery %>%
  rename(spanning.reads = num_LR)
CTATLR_Huh7_Discovery$spanning.pairs <- NA

LongGF_Huh7_Discovery <- Annot_LongGF_Huh7_Discovery %>%
  rename(spanning.reads = V3)
LongGF_Huh7_Discovery$spanning.pairs <- NA

Genion_Huh7_Discovery <- Annot_Genion_Huh7_Discovery %>%
  rename(spanning.reads = V5)
Genion_Huh7_Discovery$spanning.pairs <- NA

FusionSeeker_Huh7_Discovery <- Annot_FusionSeeker_Huh7_Discovery %>%
  rename(spanning.reads = NumSupp)
FusionSeeker_Huh7_Discovery$spanning.pairs <- NA

GFSeeker_Huh7_Discovery <- Annot_GFSeeker_Huh7_Discovery %>%
  rename(spanning.reads = "support num")
GFSeeker_Huh7_Discovery$spanning.pairs <- NA

JAFFAL_Huh7_Discovery <- Annot_JAFFAL_Huh7_Discovery %>%
  rename(spanning.reads = spanning.reads) 

JAFFAL3Gene_Huh7_Discovery <- JAFFAL3Gene_Huh7_Discovery %>%
  rename(spanning.reads = Reads) 
JAFFAL3Gene_Huh7_Discovery$spanning.pairs <- NA

cols <- c("library_type", "Platform", "full_label", 
          "Cell_Line", "RNA_sample",
          "Algorithm",  "fusionType", "Discovery",
          "Known_Gene_Fusion", "Validation", "Sample", 
          "fusion.gene.id", "fusion.ens.gene.id", 
          "spanning.reads", "spanning.pairs")

df_list <- list(
  CTATLR_Huh7_Discovery, Genion_Huh7_Discovery, 
  LongGF_Huh7_Discovery, FusionSeeker_Huh7_Discovery, 
  GFSeeker_Huh7_Discovery, JAFFAL_Huh7_Discovery, JAFFAL3Gene_Huh7_Discovery
)

fusions_rs_Huh7_discovery <- df_list %>%
  map(~ .x %>% unique() %>% select(all_of(cols))) %>%
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
