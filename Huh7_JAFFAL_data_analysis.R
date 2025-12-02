####Huh7 Sequencing Library JAFFAL results analysis####
#Author: Ryley Dorney
#Date July 2-25
#####
#Load R libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)
library(ComplexUpset)
#Load JAFFAL data into R
Nanopore_sequencing <- read.csv("/bioinformatics/ryley/Library_Benchmark/Nanopore_Huh7_fusions/jaffa_results.csv", header = TRUE) %>% 
  mutate(library_type = case_when(grepl('FBA22517|FBA22660', sample)  ~ "direct_RNA",
                                  grepl('FBA43334|FBA62655', sample)  ~ "direct_cDNA",
                                  grepl('FAZ', sample)  ~ "PCR_cDNA"),
         RNA_sample = case_when(grepl('FBA62655|FBA22660|FAZ80247', sample)  ~ "Huh7_p9_8_7_B1", 
                                grepl('FBA22517|FAZ83542|FBA43334', sample) ~ "Huh7_p9_8_7_B2"), 
         Platform = "ONT")

Illumina_sequencing_IllHB1<- read.csv("/bioinformatics/ryley/Library_Benchmark/AGRF_Sequencing_Data/Illumina_NovaSeq/fastp_trimmed/Huh7B1_JAFFA_direct_analysis/jaffa_direct_Huh7B1.csv", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",           
         RNA_sample = "Huh7_p9_8_7_B1",
        Platform = "Illumina")
Illumina_sequencing_IllHB2<- read.csv("/bioinformatics/ryley/Library_Benchmark/AGRF_Sequencing_Data/Illumina_NovaSeq/fastp_trimmed/JAFFA_direct_analysis/Huh7B2_direct/jaffa_direct_Huh7B2.csv", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA",
         RNA_sample = "Huh7_p9_8_7_B2",
         Platform = "Illumina")
Illumina_sequencing <- rbind(Illumina_sequencing_IllHB1, Illumina_sequencing_IllHB2)

PacBio_Sequencing <- read.csv("/bioinformatics/ryley/Library_Benchmark/AGRF_Sequencing_Data/PacBio_Revio_full_length_cDNA/jaffa_results.csv", header = TRUE) %>% 
  mutate(library_type = "PCR_cDNA", 
         RNA_sample = case_when(grepl('B1', sample)  ~ "Huh7_p9_8_7_B1", 
                                grepl('B2', sample)  ~  "Huh7_p9_8_7_B2"), 
         Platform = "PacBio")

Huh7_JAFFAL <- rbind(Illumina_sequencing, Nanopore_sequencing, PacBio_Sequencing)#import tri-gene results
PacBio_Sequencing_3Gene <- rbind(mutate(read.table("/bioinformatics/ryley/Library_Benchmark/AGRF_Sequencing_Data/PacBio_Revio_full_length_cDNA/PB1_FLNCcDNA.fastq/PB1_FLNCcDNA.fastq.3gene_summary", header = TRUE), sample = "PB1_FLNCcDNA.fastq") , 
                                 mutate(read.table("/bioinformatics/ryley/Library_Benchmark/AGRF_Sequencing_Data/PacBio_Revio_full_length_cDNA/PB2_FLNCcDNA.fastq/PB2_FLNC.fastq.3gene_summary", header = TRUE), sample = "PB2_FLNCcDNA.fastq" ))%>% 
  mutate(library_type = "PCR_cDNA", 
         RNA_sample = case_when(grepl('B1', sample)  ~ "Huh7_p9_8_7_B1", 
                                grepl('B2', sample)  ~  "Huh7_p9_8_7_B2"), 
         Platform = "PacBio")

##load tri gene results
Tri_gene <- list.files(path="/bioinformatics/ryley/Library_Benchmark", pattern = "3gene_summary", full.names = TRUE, recursive = TRUE) 
Nano_3gene <- Tri_gene[grepl("Nanopore_Huh7_fusionschopper_filtered", Tri_gene)]
name <- 0
for (location in Nano_3gene){
  name <- name + 1
  assign(paste0("file_", name), mutate(read.table(location, header = TRUE), sample = location, library_type = case_when(grepl('FBA22517|FBA22660', sample)  ~ "direct_RNA",
                                                                                                                        grepl('FAZ', sample)  ~ "PCR_cDNA",
                                                                                                                        grepl('FBA62655|FBA43334', sample)  ~ "direct_cDNA"),
                                       RNA_sample = case_when(grepl('FAZ80247|FBA62655|FBA22660', sample)  ~ "Huh7_p9_8_7_B1",
                                                              grepl('FBA22517|FAZ83542|FBA43334', sample) ~ "Huh7_p9_8_7_B2"),
                                       Platform = "ONT"))
} 

Nanopore_sequencing_3Gene <- rbind(file_1, file_2, file_3, file_4, file_5, file_6)

Huh7_JAFFAL_3Gene <- rbind(PacBio_Sequencing_3Gene, Nanopore_sequencing_3Gene)

Huh7_JAFFAL$library_type <- factor(
  Huh7_JAFFAL$library_type,
  levels = c("direct_RNA","direct_cDNA", "PCR_cDNA")
)
Huh7_JAFFAL$Platform <- factor(
  Huh7_JAFFAL$Platform,
  levels = c("ONT", "PacBio", "Illumina")
)
Huh7_JAFFAL$RNA_sample <- factor(
  Huh7_JAFFAL$RNA_sample,
  levels = c("Huh7_p9_8_7_B2", "Huh7_p9_8_7_B1"))

#filtering out rRNA fusions and low read support fusions
obvious_library_chimeras <- Huh7_JAFFAL %>%
  filter(!grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom1) | #non-genomic chromosomes and scaffolds
           !grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom2)) 
lowread_Huh7_JAFFAL <- Huh7_JAFFAL %>% filter((spanning.pairs==0 & spanning.reads==1) |(spanning.pairs==1 & spanning.reads==0)) #fusions supported by only one spanning read or one spanning pair


filtered_Huh7_JAFFAL <- Huh7_JAFFAL %>%
  filter(grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom1) &
           grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)", chrom2)) %>%
  filter((spanning.pairs>=0 & spanning.reads>=2) |(spanning.pairs>=2 & spanning.reads>=0)) %>% 
  mutate(fusiontype = case_when(chrom1 == chrom2 ~"Intra-chromosomal", chrom1 != chrom2 ~ "Inter-chromosomal"))

filtered_Huh7_JAFFAL_3Gene <- Huh7_JAFFAL_3Gene %>%
  filter(Reads>=2)

# Define colors for other platforms
platform_colours <- c(
  "Illumina.PCR_cDNA.B1" = "#FFA500",
  "Illumina.PCR_cDNA.B2" = "#CC8400",
  "PacBio.PCR_cDNA.B1"   = "#C80972",
  "PacBio.PCR_cDNA.B2"   = "#8B064F",
  "ONT.direct_cDNA.B1" = "#009ACD",
  "ONT.direct_cDNA.B2" = "#00688B",
  "ONT.direct_RNA.B1"  = "#00CDCD",
  "ONT.direct_RNA.B2"  = "#008B8B",
  "ONT.PCR_cDNA.B1" = "#1f77b4",
  "ONT.PCR_cDNA.B2" = "#15496e"
)
# Merge colors
color_map <- c(platform_colours)

#number of rejects
obvious_library_chimeras$fill_id <- with(
  obvious_library_chimeras,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)
lowread_Huh7_JAFFAL$fill_id <- with(
  lowread_Huh7_JAFFAL,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)

ggplot(obvious_library_chimeras, 
       aes(x = interaction(library_type, Platform),
           fill = fill_id)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "library chimeras", x = "Library")

ggplot(lowread_Huh7_JAFFAL, 
       aes(x = interaction(library_type, Platform),
           fill = fill_id)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "fusions with less than 2 supporting reads", x = "Library")

lowread_spanning_Huh7_JAFFAL <- filtered_Huh7_JAFFAL %>% filter(spanning.reads==1) #fusions supported by only one spanning read or one spanning pair
lowread_spanning_Huh7_JAFFAL$fill_id <- with(
  lowread_spanning_Huh7_JAFFAL,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)

ggplot(lowread_spanning_Huh7_JAFFAL, 
       aes(x = interaction(library_type, Platform),
           fill = fill_id)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  #scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "fusions with less than 2 spanning reads", x = "Library")

#plot the number of fusions
filtered_Huh7_JAFFAL <- filtered_Huh7_JAFFAL %>% filter(spanning.reads>=2)


ggplot(filter(filtered_Huh7_JAFFAL, spanning.reads>=2))+
  geom_point(aes(x=spanning.pairs, y=spanning.reads, colour=Platform, shape=library_type))

filtered_Huh7_JAFFAL$library_type <- factor(
  filtered_Huh7_JAFFAL$library_type,
  levels = c("direct_RNA","direct_cDNA", "PCR_cDNA")
)
filtered_Huh7_JAFFAL$Platform <- factor(
  filtered_Huh7_JAFFAL$Platform,
  levels = c("ONT", "PacBio", "Illumina")
)
filtered_Huh7_JAFFAL$RNA_sample <- factor(
  filtered_Huh7_JAFFAL$RNA_sample,
  levels = c("Huh7_p9_8_7_B1", "Huh7_p9_8_7_B2"))

# Create combined variable for fill
filtered_Huh7_JAFFAL$fill_id <- with(
  filtered_Huh7_JAFFAL,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)

# Plot
ggplot(filter(filtered_Huh7_JAFFAL, spanning.reads>=2), 
       aes(x = interaction(library_type, Platform),
           fill = fill_id)) +
  geom_bar(position = position_dodge(preserve = "single")) +
  #scale_y_break(c(4000, 1000000), scales = c(1, 4)) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "number of fusions", x = "Library")

#tri-gene fusions
# Create combined variable for fill
filtered_Huh7_JAFFAL_3Gene$fill_id <- with(
  filtered_Huh7_JAFFAL_3Gene,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)
ggplot(filtered_Huh7_JAFFAL_3Gene) +
  geom_bar(aes(x = interaction(library_type, Platform),
               fill = fill_id), position = position_dodge(preserve = "single"))+
  scale_fill_manual(values = color_map) +
  coord_flip() +
  labs(y="Number of tri-gene fusions")

#How many unique fusion genes were detected per sample, library type and platform?
#so I think this is ignoring breakpoints???
unique(filter(filtered_Huh7_JAFFAL, spanning.reads>=2)[c(2,19:23)])%>%
  distinct(fusion.genes, library_type, RNA_sample, Platform, fill_id) %>%
  count(RNA_sample, Platform, library_type, fill_id, name = "num_fusions")%>%
  ggplot(aes(x = interaction(library_type, Platform), y=num_fusions,
             fill = fill_id)) +
  geom_col(position = position_dodge(preserve = "single")) +
  #scale_y_break(c(3000, 700000), scales = c(1, 4)) +
  scale_y_log10()+
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "number of fusions", x = "Library")
  

#number of breakpoint options identified per fusion
filter(filtered_Huh7_JAFFAL, spanning.reads>=2) %>%
  count(library_type, RNA_sample, Platform, fill_id, fusion.genes) %>%  # count per fusion
  filter(n > 1) %>%                                                     # keep only repeats
  count(RNA_sample, Platform, library_type, fill_id, name = "num_fusions") %>%  # count again for plotting
  ggplot(aes(x = interaction(library_type, Platform),
             y = num_fusions,
             fill = fill_id)) +
  geom_col(position = position_dodge(preserve = "single")) +
  #scale_y_break(c(100, 100000), scales = c(1, 1)) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "number of fusions with more than one putative breakpoint", x = "Library")

#number of unique fusions identified
filter(filtered_Huh7_JAFFAL, spanning.reads>=2) %>%
  count(library_type, RNA_sample, Platform, fill_id, fusion.genes) %>%  # count per fusion
  filter(n == 1) %>%                                                     # keep only repeats
  count(RNA_sample, Platform, library_type, fill_id, name = "num_fusions") %>%  # count again for plotting
  ggplot(aes(x = interaction(library_type, Platform),
             y = num_fusions,
             fill = fill_id)) +
  geom_col(position = position_dodge(preserve = "single")) +
  #scale_y_break(c(100, 100000), scales = c(1, 1)) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "number of fusions one putative breakpoint", x = "Library")

#number of unique fusions identified
filter(filtered_Huh7_JAFFAL, spanning.reads>=2)  %>%
  distinct(library_type, RNA_sample, Platform, fill_id, fusion.genes) %>%      # keep only repeats
  count(RNA_sample, Platform, library_type, fill_id, name = "num_fusions") %>%  # count again for plotting
  ggplot(aes(x = interaction(library_type, Platform),
             y = num_fusions,
             fill = fill_id)) +
  geom_col(position = position_dodge(preserve = "single")) +
  #scale_y_break(c(100, 100000), scales = c(1, 1)) +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = color_map) +
  labs(y = "number of fusions", x = "Library")

#overlapping fusions between library types
Upset_Plot_data <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2)[c(2,19:21,23)] %>%
  distinct(fusion.genes, fill_id) %>%
  mutate(value = 1) %>%  # mark presence
  pivot_wider(
    names_from = fill_id,
    values_from = value,
    values_fill = 0
  )

ComplexUpset::upset(
  Upset_Plot_data,
  intersect = colnames(Upset_Plot_data)[-1],   # all sets except fusion.genes
  name = "Fusion overlaps"
    )
  
Upset_Illumina <- Upset_Plot_data[c(1:3)]
Upset_PacBio <- Upset_Plot_data[c(1,10,11)]
Upset_Nanopore <- Upset_Plot_data[c(1,5:9)]
Upset_Nanopore_PCR <- Upset_Plot_data[c(1,4,5)]
Upset_Nanopore_direct_cDNA <- Upset_Plot_data[c(1,8,9)]
Upset_Nanopore_direct_RNA <- Upset_Plot_data[c(1,6,7)]

Upset_B1 <- Upset_Plot_data[c(1,2,4,6,8,11)]
Upset_B2 <- Upset_Plot_data[c(1,3,5,7,9,10)]

ComplexUpset::upset(
  Upset_B2,
  intersect = colnames(Upset_B2)[-1],   # all sets except fusion.genes
  name = "Fusion overlaps")

library(reshape2)
library(pheatmap)

fusion_mat <- Upset_Plot_data %>%
  tibble::column_to_rownames("fusion.genes")

# Jaccard distances
library(vegan)
dist_mat <- vegdist(t(fusion_mat), method = "jaccard", binary = TRUE)

# similarity matrix
sim_mat <- 1 - as.matrix(dist_mat)

pheatmap(sim_mat,
         clustering_method = "average",
         display_numbers = TRUE,
         main = "Fusion similarity between libraries")


#most promiscuous gene partners
library(ggraph)
library(igraph)

fusion_gene_hub_graph <-  function(data, low_col = "#008B8B", high_col = "#00FFFF") {
  edges <- tibble(fusion = data) %>%
    separate(fusion, into = c("gene1", "gene2"), sep = ":", remove = TRUE)
  fusion_partners_graph <- graph_from_data_frame(edges, directed = FALSE)
  
  deg <- degree(fusion_partners_graph, mode = "all")
  V(fusion_partners_graph)$degree <- deg
  hub_nodes <- V(fusion_partners_graph)[degree >= 2]$name
  
  # Get names of neighbors of hub nodes
  neighbor_nodes <- unlist(
    neighborhood(fusion_partners_graph, order = 1, nodes = hub_nodes)
  )
  neighbor_nodes <- V(fusion_partners_graph)[neighbor_nodes]$name
  
  # Unique set of node names to keep
  keep_nodes <- unique(c(hub_nodes, neighbor_nodes))
  
  # Make subgraph
  filtered_graph <- induced_subgraph(fusion_partners_graph, vids = keep_nodes)
  
  # Plot
  ggraph(filtered_graph, layout = "fr") +
    geom_edge_link(alpha = 0.4, colour = "grey70") +
    geom_node_point(aes(size = degree, colour = degree)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = Inf) +
    scale_color_gradient(guide = "legend",
                         low = low_col, high = high_col) +
    theme_void()
  visIgraph(filtered_graph)
}

fusions_Huh7 <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2, library_type == "direct_RNA")$fusion.genes %>% unique()
fusion_gene_hub_graph(fusions_Huh7)

fusions_Huh7 <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2, library_type == "direct_cDNA")$fusion.genes %>% unique()
fusion_gene_hub_graph(fusions_Huh7, low_col ="#00688B", high_col = "#00BFFF")

fusions_Huh7 <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2, library_type == "PCR_cDNA", Platform == "ONT")$fusion.genes %>% unique()
fusion_gene_hub_graph(fusions_Huh7, low_col ="#15496E", high_col = "#84C3E1")

fusions_Huh7 <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2, library_type == "PCR_cDNA", Platform == "PacBio")$fusion.genes %>% unique()
fusion_gene_hub_graph(fusions_Huh7, low_col ="#8B064F", high_col = "#C80972")

fusions_Huh7 <- filter(filtered_Huh7_JAFFAL, spanning.reads>=2, library_type == "PCR_cDNA", Platform == "Illumina")$fusion.genes %>% unique()
fusion_gene_hub_graph(fusions_Huh7, low_col ="#CC8400", high_col = "#FFA500")

library(visNetwork)
visIgraph(fusion_partners_graph)


#number of previously identified fusions in Huh7

#number of previously identified fusions in HCC
#number of previously identified fusions in Liver Cancer
##Mitelman database
#number of previously identified fusions in Cancer
##Mitelman database

#number of reads for each library vs read lengths for each library vs read depth for each library (table)










