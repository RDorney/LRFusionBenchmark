##########################################
# Import Known Huh7  Fusion list
##########################################
library(readr)
KnownHuh7FusionsLiterature_2026<-read.csv("Known_Fusions/KnownHuh7FusionsLiterature_2026.csv", header = TRUE)

View(KnownHuh7FusionsLiterature_2026) 

library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("EnsDb", "Hsapiens", "110"))
ensembldbv110 <- ah[["AH113665"]]

latest_human <- query(ah, c("EnsDb", "Hsapiens"))
edb_latest <- latest_human[[max(names(latest_human))]]

standard_chrs <- c(1:22, "X", "Y", "M", "MT")

# Step 1: Map your old synonyms to IDs using the LATEST version
# (Ensembl's GENENAME often contains synonyms in its search index)
gene_name_check<- checkGeneSymbols(unique(c(KnownHuh7FusionsLiterature_2026$Gene1, KnownHuh7FusionsLiterature_2026$Gene2)),
                                   species = "human") %>%
  mutate(
    Suggested.Symbol = coalesce(Suggested.Symbol, x)
  )
modern_mapping <- ensembldb::select(edb_latest, 
                                    keys = unique(c(KnownHuh7FusionsLiterature_2026$Gene1, KnownHuh7FusionsLiterature_2026$Gene2, 
                                                    gene_name_check$Suggested.Symbol)), 
                                    keytype = "GENENAME", 
                                    columns = c("GENENAME", "GENEID", "SYMBOL", "SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))
# Step 2: Use those IDs to find what v110 calls them
v110_lookup <- ensembldb::select(ensembldbv110, 
                                 keys = modern_mapping$GENEID, 
                                 keytype = "GENEID", 
                                 columns = c("GENEID", "SYMBOL", "SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  dplyr::filter(!grepl("^LRG_", GENEID))%>%
  rename("chromosome"="SEQNAME")

dictionary <- left_join(v110_lookup,  modern_mapping, 
                        by = "GENEID") %>% 
  left_join(gene_name_check, 
            by = c("GENENAME"="Suggested.Symbol")) %>%
  dplyr::rename(original.name='x', v110_name= SYMBOL.x)%>%
  dplyr::select(c(GENEID, original.name, v110_name))

known_fusions_ah <- left_join(KnownHuh7FusionsLiterature_2026,  dictionary, 
                              by  = c("Gene1"="original.name")) %>% 
  left_join(dictionary, by  = c("Gene2"="original.name"))

#write_csv(known_fusions_ah, file="~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/KnownHuh7Fusion_ah.csv")
####################################################################################
# manually annotate where there are empty spaces by googling information
####################################################################################
huh7fusions_manualannot<- read.csv("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Known_Fusions/KnownHuh7Fusion_ah.csv")
