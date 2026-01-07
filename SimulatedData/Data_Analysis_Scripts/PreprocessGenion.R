myfiles<- list.files(path = "/bioinformatics/ryley/Algorithm_Benchmark/Adapter_porechop_trimmed/Genion", pattern = "*_genion$", full.names = TRUE)
Genion_Sim<- do.call(rbind, lapply(myfiles, function(filename) {
  if (file.info(filename)$size > 0) {
    read.table(filename) %>%
      mutate(Source = basename(filename)) %>% mutate(control = ifelse(grepl("Spiked", Source), "positive", "negative"))
  } else {
    NULL
  }
}))
og_Genion <-Genion_Sim

# Annotate
Annot_Genion_Sim <- Genion_Sim %>% separate(V1, into = c("V1.1", "V1.2", "V1.3"), "::", remove=FALSE) %>% separate(V2, into = c("V2.1", "V2.2", "V2.3"), "::", remove=FALSE) %>% separate(V8, into = c("chr1", "chr2", "chr3"), ";", remove=FALSE)
Annot_Genion_Sim[c("V1", "V2")] <- lapply(Genion_Sim[c("V1", "V2")], function(x) gsub("::", ":", x))
Annot_Genion_Sim <- Annot_Genion_Sim %>% left_join(Simulated_Fusion_Info_2, by = c('V1'='fusion.gene.id')) 
Annot_Genion_Sim$fusionType <- mapply(function(g1, g2, g3, current_type, chr1, chr2, chr3, gene_name1, gene_name2) {
  # Check if the current fusionType is empty
  if (current_type == "" || is.na(current_type)) {
    
    #check if this is two part or tri part fusion
    if(g3 == ""|| is.na(g3)){
      
      #check if part of tri fusion
      if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g2)) |
              str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g1, ":", g2)))) {
        return("truncated_tri_fusion")} 
      else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)) |
                   str_detect(paste(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V1, subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$V3, sep = ":"), paste0(g2, ":", g1)))){
        return("reverse_order:truncated_tri_fusion")} 
      
      #check for reverse order two part genes
      else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))){
        matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1)))]
        return(paste("reverse_order:", matching_row[1]))}
    } 
    #check for chromosomal misalignment
    else if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))){
      matching_row <- subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$fusionType[which(str_detect(subset(Simulated_Fusion_Info_2, fusionType != "tri_fusion")$original.fusion.gene.id, paste0(g1, ":", g2)))]
      return(paste("chromosomal_misalignment:", matching_row[1])) 
    }
    
    #check if this fusion contains mitochondrial genes 
    if ((grepl("chrM:", chr1, ignore.case = TRUE) & (!grepl("chrM:", chr2, ignore.case = TRUE) | !grepl("chrM:", chr3, ignore.case = TRUE))) |
        (grepl("chrM:", chr2, ignore.case = TRUE) & (!grepl("chrM:", chr1, ignore.case = TRUE)| !grepl("chrM:", chr3, ignore.case = TRUE)))|
        (grepl("chrM:", chr3, ignore.case = TRUE) & (!grepl("chrM:", chr1, ignore.case = TRUE)| !grepl("chrM:", chr2, ignore.case = TRUE)))){
      return("false_fusion:mitochondrial_genomic")} 
    
    else if (grepl("chrM:", chr1, ignore.case = TRUE) & grepl("chrM:", chr2, ignore.case = TRUE) & grepl("chrM:", chr3, ignore.case = TRUE)){
      return("false_fusion:mitochondrial")}
    
    #check if its a butchered tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g1, ":", g2)) |
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g3, ":", g1))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g1, ":", g3, ":", g2))| 
            str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g2, ":", g1, ":", g3)))) {
      return("wrong_order:tri_fusion")}
    
    #check if reverse order tri-fusion
    if (any(str_detect(subset(Simulated_Fusion_Info_2, fusionType == "tri_fusion")$fusion.gene.id, paste0(g3, ":", g2, ":", g1)))) {
      return("reverse_order:tri_fusion")}
    
    #if blank and none of the above, its just a false fusion
    else if (gene_name1 == gene_name2){
      return("false_fusion:self_misalignment") 
    } else {
      return("false_fusion")}
  }
  #if it's not blank, its already been labelled
  else {return(current_type)}
}, Annot_Genion_Sim$V1.1, Annot_Genion_Sim$V1.2,  Annot_Genion_Sim$V1.3, Annot_Genion_Sim$fusionType, Annot_Genion_Sim$chr1, Annot_Genion_Sim$chr2 , Annot_Genion_Sim$chr3, Annot_Genion_Sim$V2.1, Annot_Genion_Sim$V2.2)
