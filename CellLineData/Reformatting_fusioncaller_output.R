CTATLR_SGNex_Annot <- CTATLR_SGNex_msa_annot  
JAFFAL_SGNex_Annot <- JAFFAL_SGNex_msa_annot
Genion_SGNex_Annot <- Genion_SGNex_msa_annot 
FusionSeeker_SGNex_Annot <-FusionSeeker_SGNex_msa_annot 
LongGF_SGNex_Annot <-LongGF_SGNex_msa_annot 

######################
# LongGF
######################
#assign ensembl gene IDs to fusion gene partners
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = unique(c(LongGF_SGNex_Annot$Gene1, LongGF_SGNex_Annot$Gene2)),
                   mart = ensemblv110)
LongGF_SGNex_Annot <-left_join(LongGF_SGNex_Annot, gene_info, by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%
  unique()

###############################
#repeat above for Fusionseeker
###############################
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"),
                   filters = "ensembl_gene_id_version",
                   values = unique(c(FusionSeeker_SGNex_Annot$Gene1, FusionSeeker_SGNex_Annot$Gene2)),
                   mart = ensemblv110)
FusionSeeker_SGNex_Annot <- left_join(FusionSeeker_SGNex_Annot, gene_info, by  = c("Gene1"="ensembl_gene_id_version")) %>% 
  left_join(gene_info, by  = c("Gene2"="ensembl_gene_id_version")) %>%  
  unique() 

########################
#repeat above for JAFFAL
########################
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = unique(c(JAFFAL_SGNex_Annot$Gene1, JAFFAL_SGNex_Annot$Gene2)),
                   mart = ensemblv109)
JAFFAL_SGNex_Annot <-left_join(JAFFAL_SGNex_Annot, gene_info, by = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  unique() 

fsep_JAFFALcl_3Gene <- JAFFAL_SGNex_3Gene %>% filter(Reads > 1) %>% separate(Fusion, into = c('Gene1', 'Gene2', 'Gene3'), sep = ":", remove = FALSE)
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = unique(c(fsep_JAFFALcl_3Gene$Gene1, fsep_JAFFALcl_3Gene$Gene2, fsep_JAFFALcl_3Gene$Gene3)),
                   mart = ensemblv110)

fsep_JAFFALcl_3Gene <-left_join(fsep_JAFFALcl_3Gene, gene_info, by  = c("Gene1"="external_gene_name")) %>% left_join(gene_info, by  = c("Gene2"="external_gene_name"))%>% left_join(gene_info, by  = c("Gene3"="external_gene_name"))%>%
  mutate(ensembl_gene_id = coalesce(ensembl_gene_id, Gene3), ensembl_gene_id.y = coalesce(ensembl_gene_id.y, Gene2), ensembl_gene_id.x = coalesce(ensembl_gene_id.x, Gene1)) %>%  unique() 

#separate fusion genes for Genion
fsep_Gencl <- Genion_SGNex %>% separate(V1, into = c("Gene1", "Gene2", "Gene3"), "::", remove=FALSE)