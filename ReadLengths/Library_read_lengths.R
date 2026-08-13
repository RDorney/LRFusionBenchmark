#Look at read lengths for each library
Longlibrary_read_lengths <- read.delim("ReadLengths/library_read_length_summary.tsv", header = TRUE)%>%
  mutate(library_type = case_when(grepl('FBA22517|FBA22660', file)  ~ "direct_RNA", 
                                  grepl('FBA62655|FBA43334', file)  ~ "direct_cDNA", 
                                  grepl('FAZ|FLNCcDNA', file)  ~ "PCR_cDNA"), 
         RNA_sample = case_when(grepl('FAZ80247|FBA62655|FBA22660|PB1_FLNCcDNA', file)  ~ "Huh7_p9_8_7_B1",
                                grepl('FBA22517|FAZ83542|FBA62655|FBA43334|PB2_FLNCcDNA', file) ~ "Huh7_p9_8_7_B2"), 
         Platform = case_when(grepl('FBA22517|FBA22660|FAZ83542|FBA62655|FBA43334|FAZ80247', file)  ~  "ONT",
                              grepl('PB', file)  ~  "PacBio"),
         Alignment = case_when(grepl('transcriptome', file)  ~  "transcriptome", 
                               grepl('stranded', file)  ~  "stranded_genome",
                               TRUE ~ "genome"))%>%
  mutate(Library = paste(Platform, library_type, sep = "_"))%>%
  mutate(Library = factor(Library,
                               levels = c("PacBio_PCR_cDNA", "ONT_PCR_cDNA", "ONT_direct_cDNA", "ONT_direct_RNA")))

Longlibrary_read_lengths$library_type <- factor(
  Longlibrary_read_lengths$library_type,
  levels = c("direct_RNA","direct_cDNA", "PCR_cDNA")
)
Longlibrary_read_lengths$Platform <- factor(
  Longlibrary_read_lengths$Platform,
  levels = c("ONT", "PacBio", "Illumina")
)
Longlibrary_read_lengths$RNA_sample <- factor(
  Longlibrary_read_lengths$RNA_sample,
  levels = c("Huh7_p9_8_7_B1", "Huh7_p9_8_7_B2"))

Longlibrary_read_lengths$fill_id <- with(
  Longlibrary_read_lengths,
  paste(Platform, library_type, sub(".*(B[12])$", "\\1", RNA_sample), sep = ".")
)
   
Longlibrary_read_lengths$Alignment <- factor(
  Longlibrary_read_lengths$Alignment,
  levels = c("stranded_genome", "genome", "transcriptome"))

