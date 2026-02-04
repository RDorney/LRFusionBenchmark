myfiles<-list.files(path = "/bioinformatics/siyuan/gfseeker_results/simulated", 
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
  separate("break points", into = c("Breakpoint1", "Breakpoint2"), ";", remove=TRUE) %>% 
  separate(Breakpoint1, into = c("chrom1", "base1"), ":", remove=FALSE)%>%
  separate(Breakpoint2, into = c("chrom2", "base2"), ":", remove=FALSE)%>%
  mutate(Algorithm = "GFSeeker")

og_GFSeeker <- GFSeeker_SIM
