###################################
# Mito-chimeras in SGNex data
###################################
#find "Mitochondrial:Genomic" fusions chrM x genomic genes and "Mitochondrial" chimeras
mitochondrial_genomic_FusionSeeker_SGNex <- subset(FusionSeeker_SGNex, (Chrom1 == "chrM" & Chrom2 != "chrM") | (Chrom1 != "chrM" & Chrom2 == "chrM")) %>% 
  mutate(fusionType = "Mitochondrial:Genomic")
mitochondrial_genomic_LongGF_SGNex <- filter(LongGF_SGNex, (grepl("chrM:", V4, ignore.case = TRUE) & !grepl("chrM:", V5, ignore.case = TRUE)) | (!grepl("chrM:", V4, ignore.case = TRUE) & grepl("chrM:", V5, ignore.case = TRUE))) %>% 
  mutate(fusionType = "Mitochondrial:Genomic")
mitochondrial_genomic_CTATLR_SGNex <- filter(CTATLR_SGNex, (grepl("chrM:", LeftBreakpoint, ignore.case = TRUE) & !grepl("chrM:", RightBreakpoint, ignore.case = TRUE)) | (!grepl("chrM:", LeftBreakpoint, ignore.case = TRUE) & grepl("chrM:", RightBreakpoint, ignore.case = TRUE))) %>% 
  mutate(fusionType = "Mitochondrial:Genomic") #nothing observed

mito_chimeras_FusionSeeker_SGNex <- subset(FusionSeeker_SGNex, (Chrom1 == "chrM" & Chrom2 == "chrM")) %>% mutate(fusionType = "Mitochondrial")
mito_chimeras_LongGF_SGNex <- filter(LongGF_SGNex, (grepl("chrM:", V4, ignore.case = TRUE) & grepl("chrM:", V5, ignore.case = TRUE))) %>% mutate(fusionType = "Mitochondrial")
mito_chimeras_Genion_SGNex <- filter(Genion_SGNex, grepl("M:", V8, ignore.case = TRUE)) #nothing observed
mito_chimeras_JAFFAL_SGNex <- subset(JAFFAL_SGNex, (chrom1=="chrM" | chrom2=="chrM")) #nothing observed
mito_chimeras_CTATLR_SGNex <- filter(CTATLR_SGNex, (grepl("chrM:", LeftBreakpoint, ignore.case = TRUE))) %>% 
  mutate(fusionType = "Mitochondrial") #nothing observed

#############################
#combine into one data frame
#############################
mito_chimeras_SGNex <- rbind(mitochondrial_genomic_LongGF_SGNex[c(3, 6:11)], mito_chimeras_LongGF_SGNex[c(3, 6:11)])
mito_chimeras_SGNex$NumSupp <-mito_chimeras_SGNex$V3
mito_chimeras_SGNex <-mito_chimeras_SGNex[c(2:8)] %>% rbind(mito_chimeras_FusionSeeker_SGNex[c(4,10:15)], mitochondrial_genomic_FusionSeeker_SGNex[c(4,10:15)])

mito_chimeras_SGNex$Sequencing_Depth<- factor(mito_chimeras_SGNex$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
mito_chimeras_SGNex$Library<- factor(mito_chimeras_SGNex$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))
######################################
#plot the mito_chimeras_SGNex data
######################################
ggplot(filter(mito_chimeras_SGNex, NumSupp >=2), aes(x = Sequencing_Depth, fill = fusionType)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "False calls by each algorithm", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(Library+Cell_Lines~Algorithm) +
  scale_y_log10()

ggplot(filter(mito_chimeras_SGNex, NumSupp >=2), aes(x = Sequencing_Depth, y = NumSupp)) +
  geom_boxplot(aes(fill = fusionType)) +
  labs(title = "Distribution of Read Support for false chimeras",
       x = "Read Depth",
       y = "Read Support") +
  theme_minimal() +
  facet_grid(Library+Cell_Lines~ Algorithm) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####################################
# Self Aligned Fusions
#####################################
self_alignment_LongGF_SGNex <- LongGF_SGNex %>% 
  separate(V2, into = c("V2.1","V2.2"), sep=":", remove = FALSE) %>%  
  subset(V2.1 == V2.2) %>% 
  rename("NumSupp"="V3")
self_alignment_FusionSeeker_SGNex <- FusionSeeker_SGNex %>%
  subset(Gene1 == Gene2)
self_alignment_Genion_SGNex <- Genion_SGNex %>%
  separate(V1, into = c("V1.1","V1.2"), sep="::", remove = FALSE) %>%  
  subset(V1.1 == V1.2) 
self_alignment_JAFFAL_SGNex <- JAFFAL_SGNex %>% 
  separate(fusion.genes, into = c("V1.1","V1.2"), sep=":", remove = FALSE) %>% 
  subset(V1.1 == V1.2)
self_alignment_CTATLR_SGNex <- CTATLR_SGNex %>% subset(LeftGene == RightGene)
self_alignment_results <- rbind(self_alignment_LongGF_SGNex[c(5, 8:12)], 
                                self_alignment_FusionSeeker_SGNex[c(4,10:14)]) %>% 
  mutate(fusionType = "Self-Misalignment")

###############################
# All false chimeras
###############################

#combine paralogous gene fusions and mitochondrial chimeras into a dataframe containing obvious false_chimeras
false_chimeras <- rbind(self_alignment_results, mito_chimeras_cell_line, )
false_chimeras$Sequencing_Depth<- factor(false_chimeras$Sequencing_Depth, levels = c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb", "Total"))
false_chimeras$Library<- factor(false_chimeras$Library, levels = c("direct-RNA", "direct-cDNA", "PCR-cDNA"))

supp_Figure14 <-ggplot(filter(false_chimeras, NumSupp >=2, Sequencing_Depth == c("1Gb", "2.5Gb", "5Gb", "7.5Gb", "10Gb")), aes(x = Sequencing_Depth, fill = fusionType)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "False calls by each algorithm", subtitle = "minimum read support of 2", x = "Read Depth", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_grid(Library+Cell_Lines~Algorithm) +
  scale_y_log10()+
  labs(fill = "Fusion Type") 
ggsave("supp_Figure14.png", plot = supp_Figure14,  width = 210, height =297  , unit ="mm")





