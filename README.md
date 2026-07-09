# LongReadFusionCallerBenchmark

This repository contains code, workflows, and analysis files used to evaluate multiple algorithms for identifying **gene fusions in long-read RNA sequencing data**.  
Benchmarking is performed across multiple **sequencing depths** and **read quality levels**.

The project includes benchmarking of **real and simulated Oxford Nanopore Technologies (ONT) data**, with a focus on fusion detection accuracy, breakpoint precision, false-positive rates, and cross-algorithm comparisons.

---

## 🔬 Project Overview

Long-read RNA sequencing captures full-length transcripts, making it a powerful approach for detecting gene fusions and complex transcript structures.  
This repository aims to:

- Compare long-read RNA library preparation and analysis strategies  
- Benchmark fusion detection tools:
  - **JAFFAL**
  - **LongGF**
  - **FusionSeeker**
  - **Genion**
- Evaluate breakpoint accuracy and fusion classification
- Quantify false-positive and partial fusion calls
- Characterize algorithm performance across:
  - Sequencing depth (1–100 GB)
  - Read quality (e.g. Q85–Q95)
  - Simulated vs real cell line datasets
- Provide reproducible pipelines and analysis scripts for long-read fusion benchmarking

---

## 📁 Repository Structure

### 📂 `CellLineData/`

Analysis of **real cell line datasets** (K562 and MCF7).

**Key contents:**

- R Markdown workflows for algorithm benchmarking:
  - `2025_Algorithm_Cell_lines_K562_MCF7_analysis.rmd`
  - `K562_MCF7_analysis.rmd`
- Fusion preprocessing and filtering scripts:
  - `PreprocessSGNexFusions.R`
  - `CheckAtypicalFusions.R`
  - `Mito_Chimeras_SGNex.R`
  - `LongGF_extractSumGF.sh`
- Curated fusion references:
  - `known_fusions_manual_annotation.xlsx`
  - `real_gene_fusions.tsv`
- Summary figures and reports:
  - `Cell_Line_novelvsknown.pdf`

---

### 📂 `Figures/`

Final and supplementary figures used in the manuscript, along with scripts used to generate them.

- Manuscript and supplementary figures (`.pdf`, `.png`, `.jpeg`)
- Figure-generation scripts

#### 📂 `Figures/SCRIPTS/`

Modular R scripts used to generate specific plots and summaries, including:

- Fusion recall, precision, and F1 analyses
- Breakpoint statistics
- Fusion type overlap (UpSet plots)
- Read-support comparisons
- Manual filtering and true-call overlap analyses

---

### 📂 `SimulatedData/`

All resources related to **simulated fusion benchmarking**.

#### 📂 `SimulatedData/Data_Analysis_Scripts/`

Scripts for preprocessing, consolidation, and benchmarking simulated datasets:

- Preprocessing pipelines for JAFFAL, Genion, FusionSeeker, and LongGF
- Precision/recall and read-support analyses
- Breakpoint accuracy and overlap statistics
- Resource usage and performance monitoring scripts
- R Markdown analysis:
  - `simulated_reads_recall.rmd`

#### 📂 `SimulatedData/FASTQC_files/`

FastQC reports for simulated and spike-in datasets across multiple:

- Sequencing depths
- Read quality thresholds

#### 📂 `SimulatedData/SantiyCheck/`

Quality-control checks to validate spike-in recovery:

- `Check_Spiked_in_Fusions.Rmd`

#### 📂 `SimulatedData/simulated_data_info/`

Metadata and reference files describing simulated fusion datasets:

- `FUSIM_benchmark_fusions.txt`  
  *(FUSIM output describing simulated fusion events)*
- `Fusion_ReadCounts.csv`  
  *(Read counts per fusion per simulation)*
- `simulated_fusion_info_with_GeneID.tsv`
- `Spiked_Fusions.csv`
- `Spiked_Fusions_no_chimeras.csv`
- `Spiked_Fusions_with_chimeras.csv`
- `Porechop_simulate_files.rmd`

---

