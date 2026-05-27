#########################################################
# Check how many fusions were shared by multiple callers
#
#########################################################

#################
# Prep data
#################
Overlap_CellData_JAFFAL <- full_join(disc_JAFFAL_SGNex, disc_JAFFAL_SGNex_3Gene )%>%
  rename(read_support = spanning.reads) %>%
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

Overlap_CellData_Genion <- disc_Genion_SGNex %>% 
  rename(fusionGeneID = V1, read_support = V5)%>% 
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

Overlap_CellData_FusionSeeker <- disc_FusionSeeker_SGNex %>%
  rename(read_support = NumSupp)%>%
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

Overlap_CellData_LongGF <- disc_LongGF_SGNex %>%
  #mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%
  rename(read_support = V3)%>% 
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

Overlap_CellData_GFSeeker <- disc_GFSeeker_SGNex %>%
  #mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%
  rename(read_support = "support num")%>% 
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

Overlap_CellData_CTATLR <- disc_CTATLR_SGNex %>%
  mutate(fusionGeneID = paste0(ensembl_gene_id.x, "::", ensembl_gene_id.y)) %>%
  rename(read_support = num_LR)%>% 
  select(c(fusionGeneID, read_support, Source, Cell_Lines, Sequencing_Depth, Algorithm, Library))

###############################
# Combine Data
###############################
head(rbind(Overlap_CellData_JAFFAL, Overlap_CellData_Genion, 
           Overlap_CellData_FusionSeeker, Overlap_CellData_LongGF,
           Overlap_CellData_GFSeeker, Overlap_CellData_CTATLR))

fusion_data <- rbind(Overlap_CellData_JAFFAL, Overlap_CellData_Genion, 
                     Overlap_CellData_FusionSeeker, Overlap_CellData_LongGF,
                     Overlap_CellData_GFSeeker, Overlap_CellData_CTATLR)

Overlap_CellData_algorithms <- fusion_data %>% 
  group_by(fusionGeneID, Cell_Lines, Sequencing_Depth, Library) %>%
  summarise(overlap = n_distinct(Algorithm), .groups = "drop") %>%
  left_join(fusion_data, 
            by = c("fusionGeneID", "Cell_Lines", "Sequencing_Depth", "Library"))%>%
  group_by(fusionGeneID, Cell_Lines, Sequencing_Depth, Library, Algorithm, overlap) %>%
  summarise(Total_Read_Supp = sum(read_support), .groups = "drop")
Overlap_CellData_algorithms$Sequencing_Depth<- factor(Overlap_CellData_algorithms$Sequencing_Depth, 
                                                      levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
Overlap_CellData_algorithms$Library<- factor(Overlap_CellData_algorithms$Library, 
                                             levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))

Overlap_CellData_discovery <- left_join(Overlap_CellData_algorithms, known_fusions_manual_annotation, 
                                        by = c("fusionGeneID", "Cell_Lines"="Cell_Line")) %>%  
  mutate(discovery = if_else(is.na(discovery), "novel", discovery))



# Step 1: Summarize Read_Supp per algorithm
OverlapknownOnly <- Overlap_CellData_discovery %>%
  filter(discovery != "novel") %>%
  group_by(fusionGene, Cell_Lines, Sequencing_Depth, Library, Algorithm) %>%
  summarise(Read_Supp = sum(Total_Read_Supp), .groups = "drop")

# Step 2: Compute overlap (distinct algorithms per fusion call)
Overlap_CellData_info <- OverlapknownOnly %>%
  group_by(fusionGene, Cell_Lines, Sequencing_Depth, Library) %>%
  summarise(overlap = n_distinct(Algorithm), .groups = "drop")

# Step 3: Join back to retain Algorithm column
OverlapknownOnly <- OverlapknownOnly %>%
  left_join(Overlap_CellData_info, by = c("fusionGene", "Cell_Lines", "Sequencing_Depth", "Library")) %>%
  mutate(fusionGene = fct_reorder(fusionGene, overlap, .desc = TRUE))

##########################
# Plot results
##########################
ggplot(filter(OverlapknownOnly, Sequencing_Depth == "Total"), 
       aes(x = fusionGene, y = Algorithm, fill = Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Read_Supp), color = "black", size = 3, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="", fill = "Read Support") +
  theme_minimal() +
  facet_grid(Library ~ Cell_Lines, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 1000)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="grey80"))

#Supplementary figure 12####
supp_Figure12 <-  ggplot(filter(Overlap_CellData_algorithms, Sequencing_Depth == "Total", overlap > 3), 
                         aes(x = fusionGeneID, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 4, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="", fill = "") +
  theme_minimal() +
  facet_grid(Library ~ Cell_Lines, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", na.value = "grey", midpoint = 30)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.background = element_rect(color="black", fill="grey80"),
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("supp_Figure12.png", plot = supp_Figure12,  width = 297, height = 210 , unit ="mm")
ggsave("supp_Figure12.pdf", plot = supp_Figure12,  width = 297, height = 210 , unit ="mm")

Overlap_CellData_total <- Overlap_CellData_discovery %>% filter(Sequencing_Depth == "Total") %>%
  separate()

ggplot(filter(Overlap_CellData_total, str_count(fusionGeneID, "::") == 2), 
       aes(x = fusionGeneID, y = Algorithm, fill = Total_Read_Supp)) +
  geom_tile() +
  geom_text(aes(label = Total_Read_Supp), color = "black", size = 3, alpha =1) +  # Add numbers to the tiles, adjust color for visibility
  labs(x="", y="", fill = "") +
  theme_minimal() +
  facet_grid(Library ~ Cell_Lines, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "lemonchiffon", mid="orange" ,high = "orchid", 
                       na.value = "grey", , midpoint = 30)+
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(color="black", fill="grey80"))

MCF7_total_directcDNA <- Overlap_CellData_total %>% filter(Cell_Lines == "MCF7", Library == "direct-cDNA") %>% 
  select(c(discovery, Algorithm, fusionGeneID ))

shared_fusions_cell <- Overlap_CellData_algorithms %>% filter(Sequencing_Depth == "Total")%>%
  separate(fusionGeneID, into = c("Gene1","Gene2","Gene3"), sep="::", remove = FALSE)

gene_names_shared_fusions_cell <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                        filters = "ensembl_gene_id",
                                        values = unique(c(shared_fusions_cell$Gene1, shared_fusions_cell$Gene2, shared_fusions_cell$Gene3)),
                                        mart = ensemblv110) %>% unique()

shared_fusions_cell <- left_join(shared_fusions_cell, gene_names_shared_fusions_cell, by = c("Gene1" = "ensembl_gene_id")) %>%
  left_join(gene_names_shared_fusions_cell, by = c("Gene2" = "ensembl_gene_id")) %>%
  left_join(gene_names_shared_fusions_cell, by = c("Gene3" = "ensembl_gene_id")) %>%
  mutate(ensembl_gene1 = ifelse(is.na(external_gene_name.x) |nchar(external_gene_name.x) == 0, Gene1, external_gene_name.x), 
         ensembl_gene2 = ifelse(is.na(external_gene_name.y) |nchar(external_gene_name.y) == 0, Gene2, external_gene_name.y), 
         ensembl_gene3 = ifelse(is.na(external_gene_name) | nchar(external_gene_name) == 0, Gene3, external_gene_name)) %>%
  unite("fusionGene", c(ensembl_gene1, ensembl_gene2, ensembl_gene3), sep="::", remove= FALSE, na.rm = TRUE) 

shared_fusions_cell_upset <- shared_fusions_cell %>%
  distinct(Cell_Lines, Library, Algorithm, fusionGene) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Algorithm,
    values_from = present,
    values_fill = 0
  )

shared_fusions_cell_upset  %>%
  filter(Cell_Lines == "K562", Library == "direct-cDNA") %>%
  upset(
    intersect = c("FusionSeeker", "LongGF", "JAFFAL", "Genion"),
    name = "Fusions"
  )
shared_fusions_cell_upset  %>%
  filter(Cell_Lines == "K562", Library == "direct-cDNA") %>%
  upset(
    intersect = c("FusionSeeker", "LongGF", "JAFFAL", "Genion"),
    name = "Fusions"
  )
