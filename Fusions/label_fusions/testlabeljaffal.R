test_JAFFAL_Huh7 <- JAFFAL_sensecheck_annotated
idx <- which(is.na(test_JAFFAL_Huh7$fusionType) | test_JAFFAL_Huh7$fusionType == "")# Check if the current fusionType is empty

inter_idx <- idx[test_JAFFAL_Huh7$chrom1[idx] != test_JAFFAL_Huh7$chrom2[idx]] #label inter-chromosomal fusions
test_JAFFAL_Huh7$fusionType[inter_idx] <- "inter-chromosomal"
length(test_JAFFAL_Huh7$fusionType[inter_idx])
intra_idx <- idx[test_JAFFAL_Huh7$chrom1[idx] == test_JAFFAL_Huh7$chrom2[idx]]

length(test_JAFFAL_Huh7$fusionType[intra_idx])

length(unique(paste(test_JAFFAL_Huh7$Gene1[intra_idx], test_JAFFAL_Huh7$Gene2[intra_idx])))

