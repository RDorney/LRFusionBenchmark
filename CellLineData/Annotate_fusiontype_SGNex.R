###############################
# Label other types of fusions
###############################
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

#Make read-through/cis-splice check list
adjacent_genes_table <- genes %>%
  as.data.frame() %>%
  dplyr::select(
    gene_name,
    gene_id,
    seqnames,
    start,
    end,
    strand
  ) %>%
  dplyr::arrange(seqnames, start)%>%
  dplyr::group_by(seqnames) %>%
  dplyr::mutate(
    upstream_gene   = lag(gene_name),
    downstream_gene = lead(gene_name),
    upstream_strand = lag(strand),
    downstream_strand = lead(strand)
  ) %>%
  dplyr::ungroup()

adjacent_pairs <- adjacent_genes_table %>%
  dplyr::select(
    gene_name,
    seqnames,
    strand,
    upstream_gene,
    downstream_gene
  )

#######################
# Label CTAT-LR-Fusion 
#######################
Annot_CTATLR_SGNex <- CTATLR_SGNex_msa_annot

Annot_CTATLR_SGNex$fusionType <- mapply(function(current_type, chr1, chr2, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
  }
  return(current_type)
},  
Annot_CTATLR_SGNex$fusionType, 
Annot_CTATLR_SGNex$chrom1, Annot_CTATLR_SGNex$chrom2, 
Annot_CTATLR_SGNex$LeftGene, Annot_CTATLR_SGNex$RightGene)

#######################
# Label Genion 
#######################
Annot_Genion_SGNex <- Genion_SGNex_msa_annot

Annot_Genion_SGNex$fusionType <- mapply(function(g1, g2, g3, current_type, chr1, chr2, chr3, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2 & is.na(chr3)) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2 & is.na(chr3)){
      return("inter-chromosomal") 
    } else if (!is.na(g3)) {
      return("tri-fusion")
    }
  }
  return(current_type)
} , Annot_Genion_SGNex$V1.1, Annot_Genion_SGNex$V1.2,  Annot_Genion_SGNex$V1.3, 
Annot_Genion_SGNex$fusionType, 
Annot_Genion_SGNex$chr1, Annot_Genion_SGNex$chr2 , Annot_Genion_SGNex$chr3, 
Annot_Genion_SGNex$V2.1, Annot_Genion_SGNex$V2.2)

#######################
# Label LongGF 
#######################
Annot_LongGF_SGNex <- LongGF_SGNex_msa_annot
Annot_LongGF_SGNex$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
  }
  return(current_type)
} , Annot_LongGF_SGNex$Gene1, Annot_LongGF_SGNex$Gene2, 
Annot_LongGF_SGNex$fusionType, 
Annot_LongGF_SGNex$chromosome1, Annot_LongGF_SGNex$chromosome2)

#######################
# Label FusionSeeker 
#######################
Annot_FusionSeeker_SGNex <- FusionSeeker_SGNex_msa_annot
Annot_FusionSeeker_SGNex$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
  }
  return(current_type)
} , Annot_FusionSeeker_SGNex$Gene1, Annot_FusionSeeker_SGNex$Gene2, 
Annot_FusionSeeker_SGNex$fusionType, 
Annot_FusionSeeker_SGNex$Chrom1, Annot_FusionSeeker_SGNex$Chrom2)

#######################
# Label FusionSeeker 
#######################
Annot_GFSeeker_SGNex <- GFSeeker_SGNex_msa_annot
Annot_GFSeeker_SGNex$GFType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
  }
  return(current_type)
} , Annot_GFSeeker_SGNex$gene1_name, Annot_GFSeeker_SGNex$gene2_name, 
Annot_GFSeeker_SGNex$fusionType, 
Annot_GFSeeker_SGNex$chrom1, Annot_GFSeeker_SGNex$chrom2)

#######################
# Label JAFFAL 
#######################
Annot_JAFFAL_SGNex <- JAFFAL_SGNex_msa_annot
Annot_JAFFAL_SGNex$fusionType <- mapply(function(g1, g2, current_type, chr1, chr2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    #check for interchromsomal fusions
    if (chr1 == chr2) {
      return("intra-chromosomal")
      #check for reverse order fusions
    } else if (chr1 != chr2){
      return("inter-chromosomal") 
    } 
  }
  return(current_type)
} , Annot_JAFFAL_SGNex$Gene1, Annot_JAFFAL_SGNex$Gene2, 
Annot_JAFFAL_SGNex$fusionType, 
Annot_JAFFAL_SGNex$chrom1, Annot_JAFFAL_SGNex$chrom2)

#Repeat for tri-gene fusions
Annot_JAFFAL_SGNex_3Gene <- JAFFAL_SGNex_3Gene
Annot_JAFFAL_SGNex_3Gene$fusionType <- "tri-fusion"
   