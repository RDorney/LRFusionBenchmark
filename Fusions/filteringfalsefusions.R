################################################################
#filtering out rRNA fusions and low read support fusions
################################################################
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

################################################################
#Create library type, platform and RNA sample as different factors
################################################################
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
  levels = c("B2", "B1"))
