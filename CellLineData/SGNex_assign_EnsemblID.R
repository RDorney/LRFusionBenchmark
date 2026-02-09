################################################
#assign ensembl gene IDs to fusion gene partners
################################################
CTATLR_SGNex_Annot <- Annot_CTATLR_SGNex  
JAFFAL_SGNex_Annot <- Annot_JAFFAL_SGNex
JAFFAL_SGNex_Annot_3Gene <- JAFFAL_SGNex_3Gene
Genion_SGNex_Annot <- Annot_Genion_SGNex 
FusionSeeker_SGNex_Annot <-Annot_FusionSeeker_SGNex 
GFSeeker_SGNex_Annot <- Annot_GFSeeker_SGNex 
LongGF_SGNex_Annot <-Annot_LongGF_SGNex 

standard_chrs <- c(1:22, "X", "Y", "M", "MT")

######################
# LongGF
######################
#gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                   filters = "external_gene_name",
#                   values = unique(c(LongGF_SGNex_Annot$Gene1, LongGF_SGNex_Annot$Gene2)),
#                   mart = ensemblv110)

# If Biomart did not work,
# use AnnotationHub()

ah <- AnnotationHub()
query(ah, c("EnsDb", "Hsapiens", "110"))
ensembldbv110 <- ah[["AH113665"]]
gene_info <- ensembldb::select(ensembldbv110, 
                                 keys = unique(c(LongGF_SGNex_Annot$Gene1, LongGF_SGNex_Annot$Gene2)), 
                                 keytype = "SYMBOL", 
                                 columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
         "ensembl_gene_id"="GENEID")
LongGF_SGNex_Annot <-left_join(LongGF_SGNex_Annot, gene_info, by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%
  unique()

######################
# CTAT-LR
######################
#gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                   filters = "external_gene_name",
#                   values = unique(c(CTATLR_SGNex_Annot$LeftGene, CTATLR_SGNex_Annot$RightGene)),
#                   mart = ensemblv110)

# If Biomart did not work,
# use AnnotationHub()

gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(CTATLR_SGNex_Annot$LeftGene, CTATLR_SGNex_Annot$RightGene)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

CTATLR_SGNex_Annot <-left_join(CTATLR_SGNex_Annot, gene_info, by  = c("LeftGene"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("RightGene"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, RightGene), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, LeftGene)) %>%
  unique()

######################
# GFSeeker
######################
#gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                   filters = "external_gene_name",
#                   values = unique(c(GFSeeker_SGNex_Annot$gene1_name, GFSeeker_SGNex_Annot$gene2_name)),
#                   mart = ensemblv110)

# If Biomart did not work,
# use AnnotationHub()
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(GFSeeker_SGNex_Annot$gene1_name, GFSeeker_SGNex_Annot$gene2_name)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")


GFSeeker_SGNex_Annot <-left_join(GFSeeker_SGNex_Annot, gene_info, 
                               by  = c("gene1_name"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("gene2_name"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, gene2_name), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, gene1_name)) %>%
  unique()

###############################
#repeat above for Fusionseeker
###############################
#gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"),
#                   filters = "ensembl_gene_id_version",
#                   values = unique(c(FusionSeeker_SGNex_Annot$Gene1, FusionSeeker_SGNex_Annot$Gene2)),
#                   mart = ensemblv110)
#FusionSeeker_SGNex_Annot <- left_join(FusionSeeker_SGNex_Annot, gene_info, by  = c("Gene1"="ensembl_gene_id_version")) %>% 
#  left_join(gene_info, by  = c("Gene2"="ensembl_gene_id_version")) %>%  
#  unique() 
# If Biomart did not work,
# use AnnotationHub()
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(FusionSeeker_SGNex_Annot$Gene1, FusionSeeker_SGNex_Annot$Gene2)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

FusionSeeker_SGNex_Annot <- left_join(FusionSeeker_SGNex_Annot, gene_info, by  = c("Gene1"="ensembl_gene_id")) %>% 
  left_join(gene_info, by  = c("Gene2"="ensembl_gene_id")) %>%  
  unique() 

########################
#repeat above for JAFFAL
########################

#gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                   filters = "external_gene_name",
#                   values = unique(c(JAFFAL_SGNex_Annot$Gene1, JAFFAL_SGNex_Annot$Gene2, 
#                                     JAFFAL_SGNex_Annot_3Gene$Gene1, JAFFAL_SGNex_Annot_3Gene$Gene2, JAFFAL_SGNex_Annot_3Gene$Gene3)),
#                   mart = ensemblv109)


# If Biomart did not work,
# use AnnotationHub()
query(ah, c("EnsDb", "Hsapiens", "109"))
ensembldbv109 <- ah[["AH109606"]]
gene_info <- ensembldb::select(ensembldbv109, 
                               keys = unique(c(JAFFAL_SGNex_Annot$Gene1, JAFFAL_SGNex_Annot$Gene2, 
                                               JAFFAL_SGNex_Annot_3Gene$Gene1, JAFFAL_SGNex_Annot_3Gene$Gene2, JAFFAL_SGNex_Annot_3Gene$Gene3)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")

JAFFAL_SGNex_Annot <-left_join(JAFFAL_SGNex_Annot, gene_info, 
                               by = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  
  unique() 

JAFFAL_SGNex_Annot_3Gene <-left_join(JAFFAL_SGNex_Annot_3Gene, gene_info, 
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
