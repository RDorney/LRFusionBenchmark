# LongReadFusionCallerBenchmark

This repository contains code, workflows, and analysis files used to evaluate multiple algorithms developed for identifying gene fusions in long-read RNA sequencing data.
Benchmark is across multiple sequencing depths and quality levels.  
The project includes benchmarking of real and simulated ONT data, with a
focus on fusion accuracy, breakpoint precision, and algorithmic comparisons.

## 🔬 Project Overview

Long-read RNA sequencing provides full-length transcripts, making it a powerful
approach for detecting gene fusions and structural transcript variants.
This repository aims to:

- Compare performance of long-read RNA library preparation strategies  
- Benchmark fusion detection tools (JAFFAL, LongGF, FusionSeeker, Genion)
- Evaluate breakpoint accuracy and false-positive rates
- Characterize algorithm performance across:
  - Sequencing depth (1–100 GB)
  - Read quality (e.g., 85%, 95%)
  - Simulated vs. real datasets
- Provide reproducible pipelines and analysis scripts for fusion benchmarking

## 📁 Repository Structure
### Figures
- png files of Figures in our Manuscript
### Simulated Data
#### Data_Analysis_Scripts
- simulated_reads_recall.rmd
#### FASTQC_files
#### Sanity Check
- Checked_Spiked_in_Fusions.Rmd
#### simulated_data_info
- 
### CellLineData
- Analysis R Markdown of algorithm performance on K562 & MCF7 data 
  #### CellLineData/2025_Algorithm_Cell_lines_K562_MCF7_analysis.rmd
- List of Known Gene Fusions in K562 & MCF7
  #### /CellLineData/known_fusions_manual_annotation.xlsx
  #### /CellLineData/real_gene_fusions.tsv
