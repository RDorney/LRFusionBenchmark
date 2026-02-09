####Huh7 Sequencing Library Fusion Calling Results formatting####
#Author: Ryley Dorney
#Date Dec 30-25
#####
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
################################################################
#Load R libraries
################################################################
library(dplyr)
library(tidyr)

library(biomaRt)
library(AnnotationHub)
library(HGNChelper)

################################################################
# Prepare dictionary
################################################################
ah <- AnnotationHub()
query(ah, c("EnsDb", "Hsapiens", "110"))#Load Ensembl the version that bests correlate to gencode version 44
ensembldbv110 <- ah[["AH113665"]]
query(ah, c("EnsDb", "Hsapiens", "109"))#Load Ensembl the version that bests correlate to gencode version 43
ensembldbv109 <- ah[["AH109606"]]

################################################################
#Re-annotate fusions with geneIDs#
################################################################
gene_info <- ensembldb::select(ensembldbv109, 
                               keys = unique(c(Huh7_JAFFAL$Gene1, Huh7_JAFFAL$Gene2, 
                                               Huh7_JAFFAL_3Gene$Gene1, Huh7_JAFFAL_3Gene$Gene2, Huh7_JAFFAL_3Gene$Gene3)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")



