Cell_line_data_directory <- "~/LongReadFusionCallerBenchmark/CellLineData/"

known_fusions <- read.table(paste0(Cell_line_data_directory, "real_gene_fusions.tsv"), header=TRUE)
ensemblv110 <- useEnsembl("ensembl", 
                          dataset = "hsapiens_gene_ensembl", 
                          version = 110)
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = unique(c(known_fusions$Gene1, known_fusions$Gene2)),
                   mart = ensemblv110) %>% 
  unique()
known_fusions_biomart <- left_join(known_fusions,  gene_info, by  = c("Gene1"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("Gene2"="external_gene_name"))
write_csv(known_fusions_biomart, file="~/known_fusions_biomart.csv")

known_fusions_manual_annotation <- read_excel("~/LongReadFusionCallerBenchmark/CellLineData/known_fusions_manual_annotation.xlsx")

gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "ensembl_gene_id",
                   values = unique(c(known_fusions_manual_annotation$ensembl_gene_id.x, known_fusions_manual_annotation$ensembl_gene_id.y)),
                   mart = ensemblv110) %>% unique()

#make annotation of list of known fusions
known_fusions_biomart <-known_fusions_manual_annotation %>% left_join(gene_info, by  = c("ensembl_gene_id.x"="ensembl_gene_id")) %>% left_join(gene_info_2, by  = c("ensembl_gene_id.y"="ensembl_gene_id")) %>% select(c(1:6, 9, 12:13))%>% 
  mutate(external_gene_name.y = coalesce(external_gene_name.y, Gene2), external_gene_name.x = coalesce(external_gene_name.x, Gene1), fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y) , fusionGene = paste0(external_gene_name.x, "::", external_gene_name.y))
known_fusions_biomart$discovery <- "known"

#Depmap known fusions
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
