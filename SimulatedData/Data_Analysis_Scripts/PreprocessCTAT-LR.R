myfiles<-list.files(path = "/bioinformatics/siyuan/ctat_lr_results/simulated", 
                    pattern = "ctat-LR-fusion.fusion_predictions.abridged.tsv", full.names = TRUE, recursive = TRUE)
CTATLR_SIM<- bind_rows(
  lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 136) {
    read_tsv(filename,
               col_names = TRUE,
             show_col_types = FALSE) %>%
      mutate(Source = basename(dirname(filename)),
             control = ifelse(grepl("Spiked", Source), "positive", "negative"))
    } else {
    NULL
  }
}))

CTATLR_SIM<- CTATLR_SIM %>%
  filter(num_LR >= 2)

og_CTATLR <- CTATLR_SIM

###############################
# Check antisense sense fusions
###############################
#Load necessary library
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
#load relevant gtf file
gencodev44gtf <- import("/bioinformatics/ryley/Gencode44/reference_v44/gencode.v44.annotation.gtf")

########################################
# Make function to find antisense genes
########################################
genes <- gencodev44gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

get_antisense_overlaps <- function(gene, genes_gr, id_type = c("gene_name", "gene_id")) {
  
  id_type <- match.arg(id_type)
  
  if (id_type == "gene_id") {
    # strip version if present (ENSG00000xxxx.y)
    gene <- sub("\\..*$", "", gene)
    gene_ids <- sub("\\..*$", "", genes_gr$gene_id)
    target <- genes_gr[gene_ids == gene]
  } else {
    target <- genes_gr[genes_gr$gene_name == gene]
  }
  
  if (length(target) == 0)
    stop("Gene not found in annotation")
  
  hits <- findOverlaps(target, genes_gr, ignore.strand = TRUE)
  
  query_idx   <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  
  hit_df <- data.frame(
    target_gene   = as.character(target$gene_name[query_idx]),
    target_gene_id = as.character(target$gene_id[query_idx]),
    target_strand = as.character(strand(target)[query_idx]),
    gene_name     = as.character(genes_gr$gene_name[subject_idx]),
    gene_id       = as.character(genes_gr$gene_id[subject_idx]),
    gene_type     = genes_gr$gene_type[subject_idx],
    seqnames      = as.character(seqnames(genes_gr)[subject_idx]),
    start         = start(genes_gr)[subject_idx],
    end           = end(genes_gr)[subject_idx],
    strand        = as.character(strand(genes_gr)[subject_idx]),
    stringsAsFactors = FALSE
  )
  
  hit_df %>%
    dplyr::filter(
      gene_id != target_gene_id,   # remove self
      strand != target_strand      # opposite strand
    ) %>%
    dplyr::distinct(
      gene_name, gene_id, gene_type,
      seqnames, start, end, strand
    ) 
}

#get_antisense_overlaps("MALAT1", genes) #check function works
#get_antisense_overlaps("ENSG00000287557", genes) # check function works

###########################
# Check CTAT-LR for antisense genes
###########################
genes <- gencodev44gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
# get unique genes from both columns
genes_to_check <- unique(c(CTATLR_SIM$LeftGene, CTATLR_SIM$RightGene))

# run get_antisense_overlaps on each gene
CTATLR_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_name"))

# combine into a single data frame
CTATLR_antisense_df <- bind_rows(CTATLR_antisense_list, .id = "query_gene_index")

# add the original gene names
CTATLR_antisense_df$query_gene <- genes_to_check[as.integer(CTATLR_antisense_df$query_gene_index)]
CTATLR_antisense_df$query_gene_index <- NULL

antisense_pairs <- CTATLR_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

# Annotate
CTATLR_SIM_a1 <- CTATLR_SIM %>% separate(LeftBreakpoint, into = c("chrom1", "base1", "strand1"), ":", remove=FALSE) %>% separate(RightBreakpoint, into = c("chrom2", "base2", "strand2"), ":", remove=FALSE)
CTATLR_SIM_a1["#FusionName"] <- lapply(CTATLR_SIM_a1["#FusionName"], function(x) gsub("--", ":", x))
Gene_Names<-getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name"),
                       filters = "external_gene_name",
                       values = (unique(c(CTATLR_SIM_a1$LeftGene, CTATLR_SIM_a1$RightGene))),
                       mart = ensemblv110)
Gene_Names$chromosome_name <- paste0("chr", Gene_Names$chromosome_name)

CTATLR_SIM_a2 <- CTATLR_SIM_a1 %>% 
  left_join(Gene_Names, by= c("LeftGene" ='external_gene_name', 'chrom1'='chromosome_name')) %>%
  mutate(Gene1_ensembl_ID= coalesce(ensembl_gene_id, LeftGene)) %>% 
  dplyr::select(-ensembl_gene_id)  %>% 
  left_join(Gene_Names, by= c("RightGene" ='external_gene_name', 'chrom2'='chromosome_name')) %>%
  mutate(Gene2_ensembl_ID= coalesce(ensembl_gene_id, RightGene)) %>% dplyr::select(-ensembl_gene_id) %>%
  unite("fusion.gene.id" , c(Gene1_ensembl_ID, Gene2_ensembl_ID), sep = ":", remove= FALSE)
   
Annot_CTATLR_Sim <- CTATLR_SIM_a2 %>% 
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id', relationship = "many-to-many") 

Annot_CTATLR_Sim$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for truncated_tri_fusions
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
            str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
      return("truncated_tri_fusion")
      #check for reverse order fusions
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1))|
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
      return("reverse_order:truncated_tri_fusion") 
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
      return(paste("reverse_order:", matching_row[1])) 
      #check for chromosomal misalignment
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
      #check for false fusions
    } else if ((grepl("chrM:", chr1, ignore.case = TRUE) & !grepl("chrM:", chr2, ignore.case = TRUE)) | (!grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE))){
      return("false_fusion:mitochondrial_genomic") 
    }else if (grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE)){
      return("false_fusion:mitochondrial") 
    }else if ((gene_name1 == gene_name2) & (chr1 != chr2)){
      return("false_fusion:self_misalignment") 
    }else if (paste(g1, g2) %in% paste(antisense_pairs$query_gene, antisense_pairs$gene_name)
              | 
              paste(g2, g1) %in% paste(antisense_pairs$query_gene, antisense_pairs$gene_name)){
      return("false_fusion:Sense-Antisense") 
    }else {
      return("false_fusion")
    }
  }
  return(current_type)
}, Annot_CTATLR_Sim$Gene1_ensembl_ID, Annot_CTATLR_Sim$Gene2_ensembl_ID, Annot_CTATLR_Sim$fusionType, Annot_CTATLR_Sim$chrom1, Annot_CTATLR_Sim$chrom2, Annot_CTATLR_Sim$LeftGene, Annot_CTATLR_Sim$RightGene)


