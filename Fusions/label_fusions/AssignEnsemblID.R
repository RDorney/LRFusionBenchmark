############################################################
# Re-annotate fusions with geneIDs
# Author: Ryley Dorney
# Date Feb 26
#####
standard_chrs <-c(1:22, "X", "Y", "M", "MT")
################################################################
#Load R libraries
################################################################
library(dplyr)
library(tidyr)
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

######################
# LongGF
######################
# AnnotationHub()
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(LongGF_Huh7$Gene1, LongGF_Huh7$Gene2)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL", "SEQNAME")) %>%  
  dplyr::filter(SEQNAME %in% standard_chrs) %>%  
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id" = "GENEID")


LongGF_Huh7_ensemblID <-left_join(LongGF_Huh7, gene_info, by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%
  unique()

######################
# CTAT-LR
######################
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(CTATLR_Huh7$LeftGene, CTATLR_Huh7$RightGene)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

CTATLR_Huh7_ensemblID <-left_join(CTATLR_Huh7, gene_info, by  = c("LeftGene"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("RightGene"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, RightGene), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, LeftGene)) %>%
  unique()

######################
# GFSeeker
######################
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(GFSeeker_Huh7$gene1_name, GFSeeker_Huh7$gene2_name)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")


GFSeeker_Huh7_ensemblID <-left_join(GFSeeker_Huh7, gene_info, 
                                 by  = c("gene1_name"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("gene2_name"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, gene2_name), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, gene1_name)) %>%
  unique()

###############################
#repeat above for Fusionseeker
###############################
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(FusionSeeker_Huh7$Gene1, FusionSeeker_Huh7$Gene2)), 
                               keytype = "GENEID", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

FusionSeeker_Huh7_ensemblID <- left_join(FusionSeeker_Huh7, gene_info, by  = c("Gene1"="ensembl_gene_id")) %>% 
  left_join(gene_info, by  = c("Gene2"="ensembl_gene_id")) %>%  
  unique() 

########################
#repeat above for JAFFAL
########################
query(ah, c("EnsDb", "Hsapiens", "109"))
ensembldbv109 <- ah[["AH109606"]]

#gene_info <- ensembldb::select(ensembldbv109, 
#                               keys = unique(c(Huh7_JAFFAL$Gene1, Huh7_JAFFAL$Gene2)), 
#                               keytype = "SYMBOL", 
#                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
#  dplyr::filter(SEQNAME %in% standard_chrs)%>%
#  dplyr::rename("external_gene_name" = "SYMBOL", 
#                "ensembl_gene_id"="GENEID")

#Huh7_JAFFAL_ensemblID <-left_join(Huh7_JAFFAL, gene_info, 
#                               by = c("Gene1"="external_gene_name")) %>% 
#  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
#  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  
#  unique() 

#gene_info <- ensembldb::select(ensembldbv109, 
#                               keys = "NPIPA9", 
#                               keytype = "SYMBOL", 
#                               columns = c("GENEID", "SYMBOL","SEQNAME", "GENESEQSTART", "GENESEQEND"))%>%
#  dplyr::filter(SEQNAME %in% standard_chrs)%>%
#  dplyr::rename("external_gene_name" = "SYMBOL", 
#                "ensembl_gene_id"="GENEID")

gene_info <- ensembldb::select(ensembldbv109, 
                               keys = unique(c(Huh7_JAFFAL$Gene1, Huh7_JAFFAL$Gene2)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL", "SEQNAME", "GENESEQSTART", "GENESEQEND")) %>%
  dplyr::filter(SEQNAME %in% standard_chrs) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", "ensembl_gene_id" = "GENEID") %>%
  dplyr::filter(!grepl("^LRG_", ensembl_gene_id))%>%
  dplyr::mutate(SEQNAME = paste0("chr", SEQNAME))

Huh7_JAFFAL_ensemblID<- left_join(Huh7_JAFFAL, gene_info, 
                 by = c("Gene1"="external_gene_name", "chrom1"="SEQNAME")) %>% 
  mutate(ensembl_gene_id = coalesce(ensembl_gene_id, Gene1)) %>%  
  unique() %>%
  filter(base1 >= (GENESEQSTART - 5000), base1 <= (GENESEQEND + 5000)) %>%
  dplyr::rename("G1START" = "GENESEQSTART", "G1END" = "GENESEQEND")%>%
  left_join(gene_info, by  = c("Gene2"="external_gene_name", "chrom2"="SEQNAME"))%>%
  unique() %>%
  filter(base2 >= (GENESEQSTART - 5000), base2 <= (GENESEQEND + 5000)) %>%
  dplyr::rename("G2START" = "GENESEQSTART", "G2END" = "GENESEQEND")%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  
  unique() 

gene_info <- ensembldb::select(ensembldbv109, 
                              keys = unique(c(Huh7_JAFFAL_3Gene$Gene1, Huh7_JAFFAL_3Gene$Gene2, Huh7_JAFFAL_3Gene$Gene3)), 
                              keytype = "SYMBOL", 
                              columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

Huh7_JAFFAL_3Gene_ensemblID <-left_join(Huh7_JAFFAL_3Gene, gene_info, 
                                     by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  left_join(gene_info, by  = c("Gene3"="external_gene_name"))%>%
  mutate(ensembl_gene_id = coalesce(ensembl_gene_id, Gene3), 
         ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  
  unique() 

##################################
# Not necessary for Genion
##################################

##################################
#Arriba
##################################
Huh7_Arriba_EnsemblID <- Huh7_Arriba_filt%>%
  mutate(GENEID1=sub("\\..*", "", gene_id1),
         GENEID2=sub("\\..*", "", gene_id2))

######################
# STAR-Fusion
######################
Huh7_STARFusion_EnsemblID <-Huh7_STARFusion_filt%>%
  mutate(GENEID1=gsub(".*\\^|\\..*", "", LeftGene),
         GENEID2=gsub(".*\\^|\\..*", "", RightGene))
