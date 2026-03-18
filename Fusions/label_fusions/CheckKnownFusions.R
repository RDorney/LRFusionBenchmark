################################
# Title: Check Known Fusions
# Author: Ryley Dorney
# Date: Mar 2026
################################
# Load in Known Fusions List
source("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/ImportFormat_Known_Huh7.R")
View(huh7fusions_manualannot)
#huh7fusions_manualannot$Discovery <- "Known"
#######################
# Check CTAT-LR-Fusion
#######################
Annot_CTATLR_Huh7_Discovery<- bind_rows(
  # Forward match
  inner_join(Annot_CTATLR_Huh7, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "Known"),
  
  # Reverse match
  inner_join(Annot_CTATLR_Huh7, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
    mutate(Discovery = "Reverse Known") # or "Reverse Known" if you want to distinguish them
)%>%
  right_join(Annot_CTATLR_Huh7)%>%
  # If Discovery is still NA, it didn't find a match in either direction
  mutate(Discovery = if_else(is.na(Discovery), "Putative Novel", Discovery))
#######################
# Check JAFFA/L
#######################
JAFFAL_Huh7_Discovery <- bind_rows(
  # Forward match
  inner_join(Annot_JAFFAL_Huh7, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "Known"),
  
  # Reverse match
  inner_join(Annot_JAFFAL_Huh7, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.y", "ensembl_gene_id.y" = "GENEID.x")) %>%
    mutate(Discovery = "Reverse Known"))%>% 
  right_join(Annot_JAFFAL_Huh7)%>%
  # If Discovery is still NA, it didn't find a match in either direction
  mutate(Discovery = if_else(is.na(Discovery), "Putative Novel", Discovery))
#check 3 gene
JAFFAL3Gene_Huh7_Discovery <- bind_rows(
  # Forward match
  inner_join(Annot_JAFFAL_Huh7_3Gene, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.x" = "GENEID.x", "ensembl_gene_id.y" = "GENEID.y")) %>%
    mutate(Discovery = "contains Known 1.2"),
  
  # Reverse match
  inner_join(Annot_JAFFAL_Huh7_3Gene, huh7fusions_manualannot, 
            by = c("ensembl_gene_id.y" = "GENEID.x", "ensembl_gene_id" = "GENEID.y")) %>%
    mutate(Discovery = "contains Known 2.3") # or "Reverse Known" if you want to distinguish them
)%>%
  right_join(JAFFAL3Gene_Huh7_Discovery)%>%
  # If Discovery is still NA, it didn't find a match in either direction
  mutate(Discovery = if_else(is.na(Discovery), "Putative Novel", Discovery))

#######################
# Check LongGF
#######################
LongGF_Huh7_Known <- Annot_LongGF_Huh7
#######################
# Check FusionSeeker
#######################
FusionSeeker_Huh7_Known <- Annot_FusionSeeker_Huh7
#######################
# Check GFSeeker
#######################
GFSeeker_Huh7_Known <- Annot_GFSeeker_Huh7
#######################
# Check Genion
#######################
Genion_Huh7_Known <- Annot_Genion_Huh7