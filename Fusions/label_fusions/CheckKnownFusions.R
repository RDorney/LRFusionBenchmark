################################
# Title: Check Known Fusions
# Author: Ryley Dorney
# Date: Mar 2026
################################
# Load in Known Fusions List
source("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/ImportFormat_Known_Huh7.R")
View(KnownHuh7Fusions_ah)
#######################
# Check CTAT-LR-Fusion
#######################
CTATLR_Huh7_Known <- Annot_CTATLR_Huh7
KnownHuh7Fusions_ah %>%
#######################
# Check JAFFA/L
#######################
JAFFAL_Huh7_Known <- Annot_JAFFAL_Huh7
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