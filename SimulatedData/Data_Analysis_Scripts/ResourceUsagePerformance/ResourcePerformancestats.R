###############################
# Author : Ryley Dorney
# Date: 7 April 2026
# Title: Checking Resource Usage Stats
#######################################
library(readr)
RUP_table<- read_tsv("~/LongReadFusionCallerBenchmark/SimulatedData/Data_Analysis_Scripts/ResourceUsagePerformance/Benchmark_metrics_concat.tsv")
View(RUP_table)
library(dplyr)
filter(RUP_table, size_category != "100GB") %>%
  select(!c(dataset_type, input_size_gb)) %>%
  arrange(desc(peak_ram_gb)) %>%
  mutate(user_time_minutes = (user_time_seconds/60),
         wall_time_minutes = (wall_time_seconds/60))%>%
  filter(spiked == "yes")%>%
  View()
