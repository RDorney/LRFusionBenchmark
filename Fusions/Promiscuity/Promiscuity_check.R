# ==============================================================================
# PIPELINE: MULTI-GENE FUSION PROMISCUITY ANALYSIS (GENE & LIBRARY LEVEL)
# ==============================================================================

# Load relevant data
fusions_rs_Huh7_discovery<- read_tsv("~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/fusions_readsupport_Huh7_discovery.tsv.gz")

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(igraph)
library(ggplot2)

# --- STEP 1: PREPROCESS DATA AND RESOLVE MULTI-GENE CHAINS ---
# This step handles standard pairs, tri-gene, and four-gene fusions dynamically.

# 1. Assign unique event IDs and split concatenated Ensembl IDs into list elements
fusions_expanded <- fusions_rs_Huh7_discovery %>%
  mutate(event_id = row_number()) %>%
  mutate(all_genes = str_split(fusion.ens.gene.id, "::"))

# 2. Reshape into long format to discover all internal pairwise combinations
gene_partners_long <- fusions_expanded %>%
  unnest(all_genes) %>%
  rename(Gene = all_genes) %>%
  left_join(
    fusions_expanded %>% unnest(all_genes) %>% select(event_id, Partner = all_genes),
    by = "event_id",
    relationship = "many-to-many"
  ) %>%
  # Eliminate self-loops (e.g., GeneA partnering with GeneA)
  filter(Gene != Partner) %>%
  # Retain critical experimental design metadata
  select(library_type, Platform, RNA_sample, Algorithm, Gene, Partner, event_id)


# --- STEP 2: NEST DATA BY EXPERIMENTAL LAYOUT & COMPUTE METRICS ---
# Group the dataset by 4 target variables and create igraph networks.

library_networks <- gene_partners_long %>%
  select(library_type, Platform, RNA_sample, Algorithm, Gene, Partner) %>%
  distinct() %>% 
  group_by(library_type, Platform, RNA_sample, Algorithm) %>%
  nest()

# Calculate both library-level summary scores and gene-specific metrics
processed_results <- library_networks %>%
  mutate(
    # Convert edge lists directly into igraph network objects
    graph = map(data, ~ graph_from_data_frame(.x, directed = FALSE)),
    
    # Calculate global library summary stats
    library_metrics = map(graph, function(g) {
      d <- degree(g)
      # Handle power law fit safety for small networks
      fit_alpha <- tryCatch(power.law.fit(d)$alpha, error = function(e) NA)
      
      tibble(
        Total_Fusion_Genes    = vcount(g),
        Total_Unique_Pairs    = ecount(g),
        Network_Density       = edge_density(g),
        Avg_Partners_Per_Gene = mean(d),
        Network_Heterogeneity = if(mean(d) > 0) sd(d) / mean(d) else 0,
        Power_Law_Alpha       = fit_alpha
      )
    }),
    
    # Calculate gene-level network topology scores
    gene_metrics = map(graph, function(g) {
      tibble(
        Gene             = V(g)$name,
        Gene_Degree      = degree(g),
        Gene_Betweenness = round(betweenness(g), 2),
        Gene_Eigenvector = round(eigen_centrality(g)$vector, 3)
      )
    })
  )


# --- STEP 3: EXTRACT FLAT DATA TABLES FOR DOWNSTREAM ANALYSIS ---

# Table 1: Global evaluation of each library's promiscuity profile
library_promiscuity_summary <- processed_results %>%
  select(library_type, Platform, RNA_sample, Algorithm, library_metrics) %>%
  unnest(cols = c(library_metrics)) %>%
  ungroup()

# Table 2: Detailed lookup table to isolate specific problematic/promiscuous genes
gene_promiscuity_details <- processed_results %>%
  select(library_type, Platform, RNA_sample, Algorithm, gene_metrics) %>%
  unnest(cols = c(gene_metrics)) %>%
  ungroup()


# --- STEP 4: PRINT DATA AND PLOT GLOBAL COMPARISONS ---

# Display the main library summary comparison matrix in the console
print("--- Library-Wide Promiscuity Summary Profiles ---")
print(library_promiscuity_summary)

# Generate a multi-panel visual assessment layout
ggplot(library_promiscuity_summary, 
       aes(x = library_type, y = Avg_Partners_Per_Gene, fill = Platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  facet_grid(RNA_sample ~ Algorithm) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(
    title = "Global Library Promiscuity Assessment",
    subtitle = "Comparing average partner counts across experimental profiles",
    x = "Library Preparation Type",
    y = "Average Number of Partners per Gene",
    fill = "Sequencing Platform"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

ggplot(library_promiscuity_summary, 
       aes(x = library_type, y = Avg_Partners_Per_Gene, fill = Platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  facet_grid(RNA_sample ~ Algorithm) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(
    title = "Global Library Promiscuity Assessment",
    subtitle = "Comparing average partner counts across experimental profiles",
    x = "Library Preparation Type",
    y = "Average Number of Partners per Gene",
    fill = "Sequencing Platform"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )



library_promiscuity_profiles <- gene_partners_long %>%
  select(library_type, Platform, RNA_sample, Algorithm, Gene, Partner) %>%
  distinct() %>% 
  group_by(library_type, Platform, RNA_sample, Algorithm) %>%
  nest() %>%
  mutate(
    # Model library as a network where Genes = Nodes, Distinct Fusions = Edges
    graph = map(data, ~ graph_from_data_frame(.x, directed = FALSE)),
    
    # Calculate the requested network properties
    metrics = map(graph, function(g) {
      
      # Node-level fusion degree calculation
      d <- degree(g)
      total_genes <- vcount(g)
      total_fusions <- ecount(g)
      
      # Handle fallback empty states safely
      if (total_genes == 0 || total_fusions == 0) {
        return(tibble(
          Mean_Fusion_Degree = 0, Max_Fusion_Degree = 0,
          Prop_Promiscuous_GE2 = 0, Prop_Promiscuous_GE3 = 0, Prop_Promiscuous_GE5 = 0,
          Fusion_Degree_Entropy = 0,
          Hub_Dominance_Top1 = 0, Hub_Dominance_Top5 = 0, Hub_Dominance_Top10 = 0
        ))
      }
      
      # 1. Mean and 2. Maximum Fusion Degree
      mean_deg <- mean(d)
      max_deg  <- max(d)
      
      # 3. Proportion of Promiscuous Genes (Thresholds: >=2, >=3, >=5 partners)
      p_ge2 <- sum(d >= 2) / total_genes
      p_ge3 <- sum(d >= 3) / total_genes
      p_ge5 <- sum(d >= 5) / total_genes
      
      # 4. Fusion Degree Entropy (Distribution Profile)
      deg_prob <- d / sum(d)
      deg_entropy <- -sum(deg_prob * log2(deg_prob), na.rm = TRUE)
      
      # 5. Hub Dominance Index (Fraction of total edges touched by top 1, 5, 10 nodes)
      # In an undirected igraph graph, the sum of degrees of selected nodes divided by 
      # (2 * total_edges) represents the exact fraction of total fusion interactions they touch.
      sorted_degrees <- sort(d, decreasing = TRUE)
      
      hd_top1  <- min(sum(sorted_degrees[1]), 2 * total_fusions) / (2 * total_fusions)
      hd_top5  <- min(sum(sorted_degrees[1:min(5, total_genes)]), 2 * total_fusions) / (2 * total_fusions)
      hd_top10 <- min(sum(sorted_degrees[1:min(10, total_genes)]), 2 * total_fusions) / (2 * total_fusions)
      
      # Construct structured profile table row
      tibble(
        Mean_Fusion_Degree     = round(mean_deg, 3),
        Max_Fusion_Degree      = as.integer(max_deg),
        Prop_Promiscuous_GE2   = round(p_ge2, 4),
        Prop_Promiscuous_GE3   = round(p_ge3, 4),
        Prop_Promiscuous_GE5   = round(p_ge5, 4),
        Fusion_Degree_Entropy  = round(deg_entropy, 3),
        Hub_Dominance_Top1     = round(hd_top1, 4),
        Hub_Dominance_Top5     = round(hd_top5, 4),
        Hub_Dominance_Top10    = round(hd_top10, 4)
      )
    })
  ) %>%
  select(-data, -graph) %>%
  unnest(cols = c(metrics)) %>%
  ungroup()

View(library_promiscuity_profiles)
ALL_VALUES <- full_join(library_promiscuity_profiles, library_promiscuity_summary)
View(ALL_VALUES)
write_tsv(ALL_VALUES, file = "~/LibraryBenchmarkAnalysis/LibraryBenchmarkAnalysis_RProject/Fusions/Promiscuity/library_promiscuity_profiles.tsv.gz")

