###############################
# Check antisense sense fusions
###############################
#Load necessary library
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggsci)
#load relevant gtf file
gencodev44gtf <- import("/bioinformatics/ryley/Gencode44/reference_v44/gencode.v44.annotation.gtf")
# Prep files
CTATLR_mitocheck <- CTATLR_Huh7_ensemblID
JAFFAL_mitocheck <- Huh7_JAFFAL_ensemblID
Genion_mitocheck <- Genion_Huh7
FusionSeeker_mitocheck <- FusionSeeker_Huh7_ensemblID
LongGF_mitocheck <- LongGF_Huh7_ensemblID
GFSeeker_mitocheck <- GFSeeker_Huh7_ensemblID

###########################
# Check CTAT-LR
###########################
# Annotate with mitochondrial and self fusions
CTATLR_mitocheck_annotated <- CTATLR_mitocheck %>%
  dplyr::mutate(
    fusionType = dplyr::case_when(
      grepl("chrM:", LeftBreakpoint,  ignore.case = TRUE) &
        grepl("chrM:", RightBreakpoint, ignore.case = TRUE)
      ~ "Mitochondrial:Mitochondrial",
      xor(
        grepl("chrM:", LeftBreakpoint,  ignore.case = TRUE),
        grepl("chrM:", RightBreakpoint, ignore.case = TRUE)
      ) ~ "Mitochondrial:Genomic",
      (LeftGene == RightGene) ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

###########################
# Check Genion
###########################
# Annotate with mitochondrial and self fusions
Genion_mitocheck_annotated <- Genion_mitocheck %>%
  dplyr::mutate(
    fusionType = dplyr::case_when(
      grepl("M:", chr1,  ignore.case = TRUE) &
        grepl("M:", chr2, ignore.case = TRUE) 
      ~ "Mitochondrial:Mitochondrial",
      (
        grepl("M:", chr1,  ignore.case = TRUE) +
          grepl("M:", chr2, ignore.case = TRUE) +
          grepl("M:", chr3, ignore.case = TRUE)+
          grepl("M:", chr4, ignore.case = TRUE)
      ) == 1 ~ "Mitochondrial:Genomic",
      (V1.1 == V1.2) ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

###########################
# Check FusionSeeker
###########################
# Annotate with mitochondrial and self fusions
FusionSeeker_mitocheck_annotated <- FusionSeeker_mitocheck %>%
  dplyr::mutate(
    fusionType = dplyr::case_when(
      Chrom1 == "chrM" & Chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == Chrom1 , "chrM" == Chrom2) ~ "Mitochondrial:Genomic",
      Gene1 == Gene2 ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

###########################
# Check LongGF
###########################
# Annotate with mitochondrial and self fusions
LongGF_mitocheck_annotated <- LongGF_mitocheck %>%
  dplyr::mutate(
    fusionType = case_when(
      chromosome1 == "chrM" & chromosome2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chromosome1, "chrM" == chromosome2) ~ "Mitochondrial:Genomic",
      Gene1 == Gene2 ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

##################################
# Check GFSeeker
##################################
# Annotate with mitochondrial and self fusions
GFSeeker_mitocheck_annotated <- GFSeeker_mitocheck %>%
  dplyr::mutate(
    fusionType = dplyr::case_when(
      chrom1 == "chrM" & chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chrom1, "chrM" == chrom2) ~ "Mitochondrial:Genomic",
      gene1_name == gene2_name ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

##################################
# Check JAFFAL and update function
##################################
# Annotate with mitochondrial and self fusions
JAFFAL_mitocheck_annotated <- JAFFAL_mitocheck %>%
  mutate(
    fusionType = case_when(
      chrom1 == "chrM" & chrom2 == "chrM" ~ "Mitochondrial:Mitochondrial",
      xor("chrM" == chrom1, "chrM" == chrom2) ~ "Mitochondrial:Genomic",
      (Gene1 == Gene2) ~ "Self-Misalignment",
      TRUE ~ NA 
    )
  )

###########################
# combined dataframes
###########################
mitocheck_fusions <- rbind(dplyr::select(CTATLR_mitocheck_annotated, 
                                         c("RNA_sample", "Platform", "Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")),
                           dplyr::select(Genion_mitocheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")),
                           dplyr::select(LongGF_mitocheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")), 
                           dplyr::select(FusionSeeker_mitocheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")),
                           dplyr::select(GFSeeker_mitocheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")),
                           dplyr::select(JAFFAL_mitocheck_annotated, 
                                         c("RNA_sample",  "Platform","Cell_Line", 
                                           "Algorithm", "library_type", 
                                           "fusionType","full_label")))

################
# plots
################
mitocheck_fusions$Platform <- factor(mitocheck_fusions$Platform, levels = c("Illumina", "PacBio", "ONT")) 

mitocheck_fusions$library_type <- factor(mitocheck_fusions$library_type, levels = c("PCR_cDNA", "direct_cDNA", "direct_RNA")) 
#mitocheck_fusions$full_label <- interaction(mitocheck_fusions$Platform,mitocheck_fusions$library_type, mitocheck_fusions$RNA_sample)
mitocheck_fusions$full_label <- factor(mitocheck_fusions$full_label, levels = c("Illumina.PCR_cDNA.B1", "Illumina.PCR_cDNA.B2", 
                                                  "PacBio.PCR_cDNA.B1", "PacBio.PCR_cDNA.B2", 
                                                  "ONT.PCR_cDNA.B1","ONT.PCR_cDNA.B2", 
                                                  "ONT.direct_cDNA.B1", "ONT.direct_cDNA.B2", 
                                                  "ONT.direct_RNA.B1", "ONT.direct_RNA.B2"))


ggplot(dplyr::filter(mitocheck_fusions, !is.na(fusionType)), 
       aes(x = RNA_sample, colour = Algorithm)) +
  geom_point(stat = "count") +
  theme_minimal() +
  labs(title = "atypical fusions", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  labs(fill = "Algorithm")+
  scale_y_log10()+
  scale_colour_manual(values = Alg_colour_map)


ggplot(dplyr::filter(mitocheck_fusions, !is.na(fusionType)), 
       aes(x = RNA_sample, colour = full_label, shape = Algorithm)) +
  geom_point(stat = "count") +
  theme_minimal() +
  labs(title = "atypical fusions", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(fusionType~library_type+Platform) +
  #labs(fill = "Fusion Type")+
  scale_y_log10()+
  scale_colour_manual(values = platformlibsamp_colourmap)
