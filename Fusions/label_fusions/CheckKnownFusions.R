################################
# Title: Check Known Fusions
# Author: Ryley Dorney
# Date: Mar 2026
################################
# Load in Known Fusions List
source("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/ImportFormat_Known_Huh7.R")
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

fusions_Huh7_discovery <- df_list %>%
  map(~ .x %>% unique() %>% select(all_of(cols))) %>%
  bind_rows()

