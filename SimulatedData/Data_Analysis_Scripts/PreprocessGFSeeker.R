myfiles<-list.files(path = "/bioinformatics/siyuan/gfseeker_results_v2/simulated", 
                    pattern = "final_validated_fusions.csv", full.names = TRUE, recursive = TRUE)
GFSeeker_SIM<- bind_rows(
  lapply(myfiles, function(filename) {
    if (file.info(filename)$size > 118) {
      read_tsv(filename,
               col_names = TRUE,
               show_col_types = FALSE) %>%
        mutate(Source = basename(dirname(filename)),
               control = ifelse(grepl("Spiked", Source), "positive", "negative"))
    } else {
      NULL
    }
  }))

GFSeeker_SIM <- GFSeeker_SIM  %>%
  filter("support num" >= 2)%>%
  rename(gene1_name = "#gene1 name")%>%
  rename(gene2_name = "gene2 name")%>%
  separate("break points", into = c("Breakpoint1", "Breakpoint2"), ";", remove=TRUE) %>% 
  separate(Breakpoint1, into = c("chrom1", "base1"), ":", remove=FALSE)%>%
  separate(Breakpoint2, into = c("chrom2", "base2"), ":", remove=FALSE)%>%
  mutate(Algorithm = "GFSeeker") 

og_GFSeeker <- GFSeeker_SIM

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

####################################
# Check GFSeeker for false antisense genes
####################################
# get unique genes from both columns
genes_to_check <- unique(c(GFSeeker_SIM$gene1_name, GFSeeker_SIM$gene2_name))

# run get_antisense_overlaps on each gene
GFSeeker_antisense_list <- lapply(genes_to_check, 
                                function(g) get_antisense_overlaps(g, genes, id_type = "gene_name"))

# combine into a single data frame
GFSeeker_antisense_df <- bind_rows(GFSeeker_antisense_list, .id = "query_gene_index")

# add the original gene names
GFSeeker_antisense_df$query_gene <- genes_to_check[as.integer(GFSeeker_antisense_df$query_gene_index)]
GFSeeker_antisense_df$query_gene_index <- NULL

antisense_pairs <- GFSeeker_antisense_df %>%
  dplyr::select(query_gene, gene_name) %>%
  distinct()

# Annotate
gene_info <- ensembldb::select(ensembldbv110, 
                               keys = unique(c(GFSeeker_SIM$gene1_name, GFSeeker_SIM$gene2_name)), 
                               keytype = "SYMBOL", 
                               columns = c("GENEID", "SYMBOL","SEQNAME"))%>%
  dplyr::filter(SEQNAME %in% standard_chrs)%>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename("external_gene_name" = "SYMBOL", 
                "ensembl_gene_id"="GENEID")


GFSeeker_SIM_ensemblID <-left_join(GFSeeker_SIM, gene_info, 
                                    by  = c("gene1_name"="external_gene_name")) %>% 
  left_join(gene_info, by  = c("gene2_name"="external_gene_name"))%>% 
  mutate(ensembl_gene_id.y = coalesce(ensembl_gene_id.y, gene2_name), 
         ensembl_gene_id.x = coalesce(ensembl_gene_id.x, gene1_name)) %>%
  unique() %>%
  unite("fusion.gene.id" , c(ensembl_gene_id.x, ensembl_gene_id.y), sep = ":", remove= FALSE)

Annot_GFSeeker_Sim <- GFSeeker_SIM_ensemblID %>% 
  left_join(Simulated_Fusion_Info_2, by = 'fusion.gene.id', relationship = "many-to-many") 

Annot_GFSeeker_Sim$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2, gene_name1, gene_name2) {
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
      return(paste0("reverse_order:", matching_row[1])) 
      #check for chromosomal misalignment
    } else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste0("chromosomal_misalignment:", matching_row[1])) 
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
}, Annot_GFSeeker_Sim$ensembl_gene_id.x, Annot_GFSeeker_Sim$ensembl_gene_id.y, 
Annot_GFSeeker_Sim$fusionType, Annot_GFSeeker_Sim$chrom1, Annot_GFSeeker_Sim$chrom2, 
Annot_GFSeeker_Sim$gene1_name, Annot_GFSeeker_Sim$gene2_name)


