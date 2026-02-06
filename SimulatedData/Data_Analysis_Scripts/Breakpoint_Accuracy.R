###############################################
# Assess the Acuracy of predicted breakpoints #
###############################################
Simulated_BreakPoint <- Simulated_Fusion_Info_3 

################################################################
# Testing the code on one fusion gene to check expected output #
################################################################
JAFFAL_BreakPoint <- Annot_JAFFAL_Sim[c(1:25, 43, 50:52)] %>%
  distinct() %>%
  filter(control == "positive",
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment",
                fusionType), spanning.reads >= 2) 
#chrom1, base1, chrom2, base2

SB <- subset(Simulated_BreakPoint, fusion.gene.id == "ENSG00000009830:ENSG00000087258")
J <- subset(JAFFAL_BreakPoint, fusion.gene.id == "ENSG00000009830:ENSG00000087258" & depth == "100GB" & Sequence_Identity == "95%")
gene_names <- str_split("ENSG00000009830:ENSG00000087258", pattern = ":")[[1]] 
gene1 <- gene_names[1]
gene2 <- gene_names[2]

SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))

JB1 <- as.numeric(J$base1) - SB1
JB2 <- as.numeric(J$base2) - SB2

list(c(gene1, gene2) , c(SB1, SB2) , c(JB1, JB2))

nrow(distinct(Breakpoint_Accuracy))

Genion_BreakPoint <- Annot_Genion_Sim[c(1:27, 37, 44:46)] %>% 
  distinct() %>% 
  filter(control == "positive", 
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", fusionType)) %>%
  separate(V8, into = c("Gene1", "Gene2", "Gene3"), sep = ";") %>% 
  separate(Gene1, into = c("Gene1.exonStarts", "Gene1.exonEnds"), sep = "-") %>% 
  separate(Gene2, into = c("Gene2.exonStarts", "Gene2.exonEnds"), sep = "-") %>% 
  separate(Gene3, into = c("Gene3.exonStarts", "Gene3.exonEnds"), sep = "-") %>%
  separate(Gene1.exonStarts, into=c("C1", "Gene1.exonStarts")) %>% 
  separate(Gene2.exonStarts, into=c("C2", "Gene2.exonStarts")) %>% 
  separate(Gene3.exonStarts, into=c("C3", "Gene3.exonStarts")) %>% 
  mutate_at(vars(C1, C2, C3), funs(paste0("chr", .)))# "C1",  "Gene1.2", "C2", "Gene2.1", "Gene2.2", "C3", "Gene3.1" 

#################################
# Assess Breakpoint differences #
#################################

JAFFAL_BreakPoint <- Annot_JAFFAL_Sim[c(1:25, 43, 50:52)] %>%
  distinct() %>%
  filter(control == "positive",
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment",
                fusionType), spanning.reads >= 2) #chrom1, base1, chrom2, base2

LongGF_BreakPoint <- Annot_LongGF_Sim[c(1:24, 38, 45:47)] %>%
  distinct() %>%
  filter(control == "positive", 
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", 
                fusionType)) #chrom1, breakpoint1, chrom2, base2

FusionSeeker_BreakPoint <- Annot_FusionSeeker_Sim[c(1:23, 33, 40:42)] %>%
  distinct() %>%
  filter(control == "positive", 
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment", 
                fusionType)) #Chrom1, Breakpoint1, Chrom2, Base2

CTAT_LR_Breakpoint <- Annot_CTATLR_Sim[c(1:25, 43, 50:52)]%>%
  distinct() %>%
  filter(control == "positive",
         !grepl("false|wrong|truncated|reverse|chromosomal_misalignment",
                fusionType), num_LR >= 2) #chrom1, base1, chrom2, base2

JAFFAL_BreakPoint_Accuracy <-  data.frame(matrix(nrow = 0, ncol = length(c('seq_depth', 'seq_id', 'gene_fusion', 'JB1', 'JB2'))))
colnames(JAFFAL_BreakPoint_Accuracy) <- c('seq_depth', 'seq_id', 'gene_fusion','JB1', 'JB2')
for (seq_depth in c('1GB' , '10GB' , '100GB')) {
  for (seq_id in c('85%' , '90%' , '95%')) {   
    for(gene_fusion in unique(c(JAFFAL_BreakPoint$fusion.gene.id))){
      SB <- subset(Simulated_BreakPoint, fusion.gene.id == gene_fusion)
      
      gene_names <- str_split(gene_fusion, pattern = ":")[[1]] 
      gene1 <- gene_names[1]
      gene2 <- gene_names[2]
      
      SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
      SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))
      
      gene3 <- gene_names[3]
      SB3 <- NA
      SB4 <- NA
      if(is.na(gene3) == FALSE){
        SB3 <- as.integer(str_split(subset(SB, alternative_ID == gene2)$exonEnds, ",")[[1]])
        SB4 <- as.integer(str_split(subset(SB, alternative_ID == gene3)$exonStarts, ",")[[1]])
      }
      
      J  <- subset(JAFFAL_BreakPoint, fusion.gene.id == gene_fusion & depth == seq_depth & Sequence_Identity == seq_id) 
      if(nrow(J) > 0){
        JB1 <- as.numeric(J$base1) - SB1
        JB2 <- as.numeric(J$base2) - SB2
      }else{
        JB1 <- "NA"
        JB2 <- "NA"
      }
      new_row <-  data.frame(seq_depth, seq_id, gene_fusion, JB1, JB2)
      colnames(new_row) <- colnames(JAFFAL_BreakPoint_Accuracy)
      JAFFAL_BreakPoint_Accuracy <- rbind(JAFFAL_BreakPoint_Accuracy, new_row)
    }
  }
}

Genion_Breakpoint_Accuracy <- data.frame(matrix(nrow = 0, ncol = length(c('seq_depth', 'seq_id', 'gene_fusion', 'GB1', 'GB2', 'GB3', 'GB4'))))
colnames(Genion_Breakpoint_Accuracy) <- c('seq_depth', 'seq_id', 'gene_fusion', 'GB1', 'GB2', 'GB3', 'GB4')
for (seq_depth in c('1GB' , '10GB' , '100GB')) {
  for (seq_id in c('85%' , '90%' , '95%')) {   
    for(gene_fusion in unique(c(Genion_BreakPoint$V1))){
      SB <- subset(Simulated_BreakPoint, fusion.gene.id == gene_fusion)
      gene_names <- str_split(gene_fusion, pattern = ":")[[1]] 
      gene1 <- gene_names[1]
      gene2 <- gene_names[2]
      
      SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
      SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))
      
      gene3 <- gene_names[3]
      SB3 <- NA
      SB4 <- NA
      if(is.na(gene3) == FALSE){
        SB3 <- as.integer(str_split(subset(SB, alternative_ID == gene2)$exonEnds, ",")[[1]])
        SB4 <- as.integer(str_split(subset(SB, alternative_ID == gene3)$exonStarts, ",")[[1]])
      } 
      G <- subset(Genion_BreakPoint,  V1 == gene_fusion & depth == seq_depth & Sequence_Identity == seq_id)
      if(nrow(G) > 0){
        
        GB1 <- as.numeric(G$Gene1.exonEnds) - SB1
        GB2 <- as.numeric(G$Gene2.exonStarts) - SB2
        GB3 <- "NA"
        GB4 <- "NA"
        if(is.na(gene3) == FALSE){
          GB3 <- toString(apply(outer(as.numeric(G$Gene2.exonEnds) , SB3, "-"), 1, function(x) x[which.min(abs(x))]))
          GB4 <- toString(apply(outer(as.numeric(G$Gene3.exonStarts) , SB4, "-"), 1, function(x) x[which.min(abs(x))]))
        }
      }else{
        GB1 <- "NA"
        GB2 <- "NA"
        GB3 <- "NA"
        GB4 <- "NA"
      }
      new_row <-  data.frame(seq_depth, seq_id, gene_fusion, GB1, GB2, GB3, GB4)
      colnames(new_row) <- colnames(Genion_Breakpoint_Accuracy)
      Genion_Breakpoint_Accuracy <- rbind(Genion_Breakpoint_Accuracy, new_row)
    }
  } 
}

LongGF_Breakpoint_Accuracy <- data.frame(matrix(nrow = 0, ncol = length(c('seq_depth', 'seq_id', 'gene_fusion', 'LB1', 'LB2'))))
colnames(LongGF_Breakpoint_Accuracy) <- c('seq_depth', 'seq_id', 'gene_fusion', 'LB1', 'LB2')
for (seq_depth in c('1GB' , '10GB' , '100GB')) {
  for (seq_id in c('85%' , '90%' , '95%')) {   
    for(gene_fusion in unique(c(LongGF_BreakPoint$fusion.gene.id))){
      SB <- subset(Simulated_BreakPoint, fusion.gene.id == gene_fusion)
      gene_names <- str_split(gene_fusion, pattern = ":")[[1]] 
      gene1 <- gene_names[1]
      gene2 <- gene_names[2]
      
      SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
      SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))
      
      gene3 <- gene_names[3]
      SB3 <- NA
      SB4 <- NA
      if(is.na(gene3) == FALSE){
        SB3 <- as.integer(str_split(subset(SB, alternative_ID == gene2)$exonEnds, ",")[[1]])
        SB4 <- as.integer(str_split(subset(SB, alternative_ID == gene3)$exonStarts, ",")[[1]])
      } 
      L <- subset(LongGF_BreakPoint, fusion.gene.id== gene_fusion & depth == seq_depth & Sequence_Identity == seq_id)
      if(nrow(L) > 0){
        LB1 <- as.numeric(L$breakpoint1) - SB1
        LB2 <- as.numeric(L$breakpoint2) - SB2
      }else{
        LB1 <- "NA"
        LB2 <- "NA"
      }
      new_row <-  data.frame(seq_depth, seq_id, gene_fusion, LB1, LB2)
      colnames(new_row) <- colnames(LongGF_Breakpoint_Accuracy)
      LongGF_Breakpoint_Accuracy <- rbind(LongGF_Breakpoint_Accuracy, new_row)
    }
  } 
}

FusionSeeker_Breakpoint_Accuracy <- data.frame(matrix(nrow = 0, ncol = length(c('seq_depth', 'seq_id', 'gene_fusion', 'FSB1', 'FSB2'))))
colnames(FusionSeeker_Breakpoint_Accuracy) <- c('seq_depth', 'seq_id', 'gene_fusion', 'FSB1', 'FSB2')
for (seq_depth in c('1GB' , '10GB' , '100GB')) {
  for (seq_id in c('85%' , '90%' , '95%')) {   
    for(gene_fusion in unique(c(FusionSeeker_BreakPoint$fusion.gene.id))){
      SB <- subset(Simulated_BreakPoint, fusion.gene.id == gene_fusion)
      FS <- subset(FusionSeeker_BreakPoint, fusion.gene.id == gene_fusion & depth == seq_depth & Sequence_Identity == seq_id)
      gene_names <- str_split(gene_fusion, pattern = ":")[[1]] 
      gene1 <- gene_names[1]
      gene2 <- gene_names[2]
      
      SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
      SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))
      
      gene3 <- gene_names[3]
      SB3 <- NA
      SB4 <- NA
      if(is.na(gene3) == FALSE){
        SB3 <- as.integer(str_split(subset(SB, alternative_ID == gene2)$exonEnds, ",")[[1]])
        SB4 <- as.integer(str_split(subset(SB, alternative_ID == gene3)$exonStarts, ",")[[1]])
      }
      if(nrow(FS) > 0){
        FSB1 <- as.numeric(FS$Breakpoint1) - SB1
        FSB2 <- as.numeric(FS$Breakpoint2) - SB2
      }else{
        FSB1 <- "NA"
        FSB2 <- "NA"
      }
      new_row <-  data.frame(seq_depth, seq_id, gene_fusion, FSB1, FSB2)
      colnames(new_row) <- colnames(FusionSeeker_Breakpoint_Accuracy)
      FusionSeeker_Breakpoint_Accuracy <- rbind(FusionSeeker_Breakpoint_Accuracy, new_row)
    }
  } 
}

CTAT_LR_Breakpoint_Accuracy <- data.frame(matrix(nrow = 0, ncol = length(c('seq_depth', 'seq_id', 'gene_fusion', 'CTATLRB1', 'CTATLRB2'))))
colnames(CTAT_LR_Breakpoint_Accuracy) <- c('seq_depth', 'seq_id', 'gene_fusion', 'CTATLRB1', 'CTATLRB2')
for (seq_depth in c('1GB' , '10GB' , '100GB')) {
  for (seq_id in c('85%' , '90%' , '95%')) {   
    for(gene_fusion in unique(c(CTAT_LR_Breakpoint$fusion.gene.id))){
      SB <- subset(Simulated_BreakPoint, fusion.gene.id == gene_fusion)
      CTATLR <- subset(CTAT_LR_Breakpoint, 
                       fusion.gene.id == gene_fusion & depth == seq_depth & Sequence_Identity == seq_id)
      gene_names <- str_split(gene_fusion, pattern = ":")[[1]] 
      gene1 <- gene_names[1]
      gene2 <- gene_names[2]
      
      SB1 <- as.integer(tail(str_split(subset(SB, alternative_ID == gene1)$exonEnds, ",")[[1]], 1))
      SB2 <- as.integer(head(str_split(subset(SB, alternative_ID == gene2)$exonStarts, ",")[[1]],1))
      
      gene3 <- gene_names[3]
      SB3 <- NA
      SB4 <- NA
      if(is.na(gene3) == FALSE){
        SB3 <- as.integer(str_split(subset(SB, alternative_ID == gene2)$exonEnds, ",")[[1]])
        SB4 <- as.integer(str_split(subset(SB, alternative_ID == gene3)$exonStarts, ",")[[1]])
      }
      if(nrow(CTATLR) > 0){
        CTATLRB1 <- as.numeric(CTATLR$base1) - SB1
        CTATLRB2 <- as.numeric(CTATLR$base2) - SB2
      }else{
        CTATLRB1 <- "NA"
        CTATLRB2 <- "NA"
      }
      new_row <-  data.frame(seq_depth, seq_id, gene_fusion, CTATLRB1, CTATLRB2)
      colnames(new_row) <- colnames(CTAT_LR_Breakpoint_Accuracy)
      CTAT_LR_Breakpoint_Accuracy <- rbind(CTAT_LR_Breakpoint_Accuracy, new_row)
    }
  } 
}

Breakpoint_Accuracy <- full_join(FusionSeeker_Breakpoint_Accuracy, LongGF_Breakpoint_Accuracy) %>%
  full_join(Genion_Breakpoint_Accuracy) %>%
  full_join(JAFFAL_BreakPoint_Accuracy) %>%
  full_join(CTAT_LR_Breakpoint_Accuracy)

Breakpoint_Accuracy[is.na(Breakpoint_Accuracy)] <- "NA"

Breakpoint_Accuracy <- Breakpoint_Accuracy %>% filter(!if_all(FSB1:JB2, ~ .x == "NA")) %>%
  pivot_longer(cols = 4:ncol(Breakpoint_Accuracy), names_to =  "Tool_Breakpoint", values_to = "Distance_from_Simulated_Breakpoint") %>%
  mutate(Algorithm = case_when(
    grepl("^GB", Tool_Breakpoint) ~ "Genion",
    grepl("^JB", Tool_Breakpoint) ~ "JAFFAL",
    grepl("^LB", Tool_Breakpoint) ~ "LongGF",
    grepl("^CTATLR", Tool_Breakpoint) ~ "CTAT-LR-Fusion",
    grepl("^FSB", Tool_Breakpoint) ~ "FusionSeeker"))%>%
  mutate(values_vector = strsplit(as.character(Distance_from_Simulated_Breakpoint), ","))%>%
  unnest(values_vector) %>%
  mutate(values_vector = as.numeric(values_vector)) %>%
  mutate(Breakpoint_number = case_when(
    grepl("B1", Tool_Breakpoint) ~ "Breakpoint_1",
    grepl("B2", Tool_Breakpoint) ~ "Breakpoint_2",
    grepl("B3", Tool_Breakpoint) ~ "Breakpoint_3",
    grepl("B4", Tool_Breakpoint) ~ "Breakpoint_4"))

Breakpoint_Accuracy <-filter(Breakpoint_Accuracy, Distance_from_Simulated_Breakpoint != "NA") %>% distinct()

#some breakpoints are before the value and have a negative value, this we transform the values to absolute distances
absolute_Breakpoint_Accuracy <- Breakpoint_Accuracy %>%
  mutate(values_vector = abs(values_vector))
write_tsv(absolute_Breakpoint_Accuracy, "~/LongReadFusionCallerBenchmark/Figures/Input_dataframes/absolute_Breakpoint_Accuracy.tsv")

