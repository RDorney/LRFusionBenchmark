Cell_line_data_directory <- "~/LongReadFusionCallerBenchmark/CellLineData/"

known_fusions <- read.table(paste0(Cell_line_data_directory, "real_gene_fusions.tsv"), header=TRUE)

################################################################
# Load Biomart
#
# if biomart doesn't work, proceed below to use annotation hub
################################################################
options(timeout = 300)  # 5 minutes
ensemblv110 <- useEnsembl("ensembl", 
                          dataset = "hsapiens_gene_ensembl", 
                          host = "https://jul2023.archive.ensembl.org")
gene_info <- distinct(rbind(getBM(attributes = c("external_gene_name", "external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = unique(c(known_fusions$Gene1, known_fusions$Gene2)),
                   mart = ensemblv110), 
                   getBM(attributes = c("external_gene_name", "ensembl_gene_id", "external_synonym"),
                    filters = "external_synonym",
                    values = unique(c(known_fusions$Gene1, known_fusions$Gene2)),
                    mart = ensemblv110)))
known_fusions_biomart <- left_join(known_fusions,  gene_info, 
                                   by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))
write_tsv(known_fusions_biomart, file="~/LongReadFusionCallerBenchmark/CellLineData/known_fusions_biomart.tsv")

###########################################################
# If Biomart did not work,
# use AnnotationHub()
###########################################################
ah <- AnnotationHub()
query(ah, c("EnsDb", "Hsapiens", "110"))
ensembldbv110 <- ah[["AH113665"]]
query(ah, c("EnsDb", "Hsapiens", "109"))
ensmebldbv109 <- ah[["AH109606"]]

latest_human <- query(ah, c("EnsDb", "Hsapiens"))
edb_latest <- latest_human[[max(names(latest_human))]]

# Step 1: Map your old synonyms to IDs using the LATEST version
# (Ensembl's GENENAME often contains synonyms in its search index)
gene_name_check<- checkGeneSymbols(unique(c(known_fusions$Gene1, known_fusions$Gene2)),
                                   species = "human") %>%
  mutate(
    Suggested.Symbol = coalesce(Suggested.Symbol, x)
  )
modern_mapping <- ensembldb::select(edb_latest, 
                                    keys = unique(c(known_fusions$Gene1, known_fusions$Gene2, gene_name_check$Suggested.Symbol)), 
                                    keytype = "GENENAME", 
                                    columns = c("GENENAME", "GENEID", "SYMBOL"))%>%
  distinct(SYMBOL, .keep_all = TRUE)
# Step 2: Use those IDs to find what v110 calls them
v110_lookup <- ensembldb::select(ensembldbv110, 
                                 keys = modern_mapping$GENEID, 
                                 keytype = "GENEID", 
                                 columns = c("GENEID", "SYMBOL"))%>%
  distinct(SYMBOL, .keep_all = TRUE)


dictionary <- left_join(v110_lookup,  modern_mapping, 
          by = "GENEID") %>% 
  left_join(gene_name_check, 
            by = c("GENENAME"="Suggested.Symbol")) %>%
  rename("original.name"= "x", "v110_name"= "SYMBOL.x")%>%
  dplyr::select(c(GENEID, original.name, v110_name))

known_fusions_ah <- left_join(known_fusions,  dictionary, 
                                   by  = c("Gene1"="original.name")) %>% 
  left_join(dictionary, by  = c("Gene2"="original.name"))
write_tsv(known_fusions_ah, file="~/LongReadFusionCallerBenchmark/CellLineData/known_fusions_ah.tsv")
####################################################################################
# manually annotate where there are empty spaces by looking at ensembl information
####################################################################################
known_fusions_manual_annotation <- read_tsv("~/LongReadFusionCallerBenchmark/CellLineData/known_fusions_manual_annotation.tsv")

#make annotation of list of known fusions
known_fusions_manual_annotation <-known_fusions_manual_annotation %>% 
  rename("ensembl_gene_id_1"="GENEID.x", "ensembl_gene_id_2"="GENEID.y") %>% 
  dplyr::select(c(1:2, 6:9)) %>%
 mutate(fusionGene = paste0(v110_name.x, "::", v110_name.y),
        fusion.gene.id = paste0(ensembl_gene_id_1, "::",ensembl_gene_id_2))
known_fusions_manual_annotation$discovery <- "known"

#####################################
#Depmap known fusions
######################################
DepMap_knownfusions <- read.csv("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/DepMap_OmicsFusionFiltered.csv", 
                                header = TRUE)
Model_Depmap <- (read.csv("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Model.csv", 
                          header = TRUE))
DepMap_knownfusions <- left_join(DepMap_knownfusions, Model_Depmap) %>% select(Gene1, Gene2, StrippedCellLineName) %>% filter((StrippedCellLineName == "K562"| StrippedCellLineName =="MCF7"))

DepMap_K562 <- filter(DepMap_knownfusions, StrippedCellLineName == "K562")%>%
  extract(Gene1, into = c("Gene1_name", "Gene1_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  extract(Gene2, into = c("Gene2_name", "Gene2_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  mutate(
    Gene1_ensembl = na_if(Gene1_ensembl, "."),
    Gene2_ensembl = na_if(Gene2_ensembl, "."))

DepMap_MCF7 <- filter(DepMap_knownfusions, StrippedCellLineName == "MCF7")%>%
  extract(Gene1, into = c("Gene1_name", "Gene1_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  extract(Gene2, into = c("Gene2_name", "Gene2_ensembl"), 
          regex = "([^ ]+) \\(([^)]+)\\)", remove = FALSE) %>%
  mutate(
    Gene1_ensembl = na_if(Gene1_ensembl, "."),
    Gene2_ensembl = na_if(Gene2_ensembl, "."))

ensembl_current <-  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
GeneIDs<-unique(rbind(getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "external_synonym", "hgnc_symbol"),
                            filters = c("external_gene_name"),
                            values = unique(c(DepMap_MCF7$Gene1_name, DepMap_MCF7$Gene2_name, DepMap_K562$Gene1_name, DepMap_K562$Gene2_name)),
                            mart = ensembl_current), 
                      getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "external_synonym", "hgnc_symbol"),
                            filters = c("external_synonym"),
                            values = unique(c(DepMap_MCF7$Gene1_name, DepMap_MCF7$Gene2_name, DepMap_K562$Gene1_name, DepMap_K562$Gene2_name)),
                            mart = ensembl_current), 
                      getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "external_synonym", "hgnc_symbol"),
                            filters = c("hgnc_symbol"),
                            values = unique(c(DepMap_MCF7$Gene1_name, DepMap_MCF7$Gene2_name, DepMap_K562$Gene1_name, DepMap_K562$Gene2_name)),
                            mart = ensembl_current)))

GeneIDs_long <- GeneIDs %>%
  pivot_longer(cols = c(external_gene_name, external_synonym), 
               names_to = "name_type", values_to = "gene_symbol") %>%
  distinct()

DepMap_K562_withID <- DepMap_K562 %>%
  left_join(GeneIDs_long, by = c("Gene1_name" = "gene_symbol")) %>%
  left_join(GeneIDs_long, by = c("Gene2_name" = "gene_symbol"), suffix = c("_Gene1", "_Gene2"))

DepMap_MCF7_withID <- DepMap_MCF7 %>%
  left_join(GeneIDs_long, by = c("Gene1_name" = "gene_symbol")) %>%
  left_join(GeneIDs_long, by = c("Gene2_name" = "gene_symbol"), suffix = c("_Gene1", "_Gene2"))

DepMap_specificCells <- rbind(DepMap_K562_withID, DepMap_MCF7_withID)
DepMap_specificCells$discovery <- "known"
