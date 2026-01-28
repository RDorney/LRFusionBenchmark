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

og_GFSeeker <- GFSeeker_SIM
