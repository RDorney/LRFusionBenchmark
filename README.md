# LRS-FusionBench

**LRS-FusionBench** is a collection of scripts, datasets, and analysis workflows used to benchmark long-read RNA sequencing fusion detection tools. The repository contains code for preprocessing sequencing data, standardising fusion caller outputs, benchmarking discovery performance, evaluating known fusions, analysing computational resource usage, and generating publication-quality figures.

This repository accompanies ongoing work benchmarking long-read fusion detection algorithms using both simulated datasets and real sequencing data, including Huh7 hepatocellular carcinoma (HCC) cell line datasets.

---

## Repository structure

```
LRS-FusionBench/
│
├── CellLineData/
│   Analysis of SGNex K562 and MCF7 cell line fusion datasets
│
├── Figures/
│   Scripts and input data used to generate publication figures
│
├── Fusions/ 
│   Huh7 Fusion preprocessing, discovery benchmarking, resource analysis,
│   fusion labelling and promiscuity analyses
│
├── Known_Fusions/
│   Compilation of known cancer fusion databases and Huh7 reference fusions
│
├── PreprocessingReads/
│   Shell scripts for adapter trimming, read filtering and preprocessing
│
├── ReadLengths/
│   Read length analyses and summary statistics
│
├── SimulatedData/
│   Scripts for analysing simulated benchmarking datasets
│
├── Huh7_JAFFAL_data_analysis.R
├── Preparing_data.rmd
└── README.md
```

---

## Repository contents

### PreprocessingReads/

Scripts for preparing sequencing reads prior to fusion detection.

* `DoradoTrimAdaptor.sh` — Adapter trimming using Dorado.
* `Fastp_trim_and_filter.sh` — Read trimming and quality filtering with Fastp.
* `Filter_NanoporeReads_Chopper.sh` — Read filtering using Chopper.
* `sequencing_data_shortcuts.sh` — Convenience commands for managing sequencing datasets.

---

### Fusions/

Contains the primary benchmarking workflows performed on our Huh7 Data.

Major components include:

* Formatting and standardisation of fusion caller outputs
* Fusion discovery benchmarking
* Fusion labelling and annotation
* Resource usage benchmarking
* Fusion promiscuity analyses
* Tool execution scripts

Key files:

* `FormattingandPreprocessFusions.R`
* `fusions_Huh7_discovery.tsv.gz`
* `fusions_readsupport_Huh7_discovery.tsv.gz`

Subdirectories include:

* `Discovery/`
* `label_fusions/`
* `Promiscuity/`
* `resource_benchmark/`
* `runtools/`

---

### Known_Fusions/

Compilation of publicly available fusion databases together with manually curated Huh7 fusion references.

Included resources:

* COSMIC Fusion
* ChimerDB
* ChimerSeq
* ChimerPub
* Mitelman Database
* DepMap Fusion data
* TCGA fusion datasets

Scripts are provided to:

* import and standardise database formats
* merge multiple databases
* generate reference fusion sets for benchmarking

---

### CellLineData/

Analysis scripts for publicly available long-read sequencing datasets.

Includes analyses of:

* SGNex datasets
* -  K562
* - MCF7

Example analyses include:

* fusion overlap
* fusion discovery and recall
* atypical fusion detection
* mitochondrial chimeras
* annotation using Ensembl identifiers

---

### SimulatedData/

Benchmarking analyses using simulated long-read RNA sequencing datasets.

Contains scripts for:

* simulated dataset analysis
* FASTQC summaries
* sanity checks
* simulation metadata

---

### ReadLengths/

Scripts analysing sequencing read length distributions.

---

### Figures/

Contains scripts and intermediate data used to generate publication figures.

---

## Workflow overview

A typical benchmarking workflow consists of:

1. Preprocess sequencing reads.
2. Run one or more long-read fusion callers.
3. Standardise fusion outputs.
4. Merge results across tools.
5. Annotate fusion characteristics.
6. Compare against known fusion databases.
7. Calculate discovery, recall and precision metrics.
8. Generate summary tables and figures.

---

## Software

Most analyses are performed in **R**.

Common R packages include:

* dplyr
* tidyr
* ggplot2
* data.table
* ComplexUpset
* ensembldb
* AnnotationHub
* BiocGenerics
* GenomicRanges

Shell scripts additionally make use of common bioinformatics tools including:

* Dorado
* Fastp
* Chopper
* samtools

---

## Data

The repository contains:

* processed benchmarking datasets
* curated known fusion reference sets
* Huh7 discovery fusion datasets
* read length summaries
* scripts for reproducing analyses

Several reference databases (e.g. COSMIC, ChimerDB and DepMap) are redistributed only for convenience. Users should ensure they comply with the respective database licensing terms.

---

## Citation

If you use this repository in your work, please cite the accompanying publication (to be added).

---

## Authors

**Ryley Dorney**
**Siyuan (Thaddeus) Wu**
**Ulf Schmitz**

ComBioLab @ James Cook University

Research focus:

* Long-read RNA sequencing
* Gene fusion detection
* Bioinformatics benchmarking
---

