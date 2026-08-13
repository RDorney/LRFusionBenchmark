# LRS-FusionBench

**LRS-FusionBench** is a collection of scripts, datasets, and analysis workflows used to benchmark long-read RNA sequencing fusion detection tools. The repository contains code for preprocessing sequencing data, standardising fusion caller outputs, benchmarking discovery performance, evaluating known fusions, analysing computational resource usage, and generating publication-quality figures.

This repository accompanies ongoing work benchmarking long-read fusion detection algorithms using both simulated datasets and real sequencing data, including Huh7 hepatocellular carcinoma (HCC) cell line datasets.

## Scope

This repository documents the analysis code behind the study. The R analysis scripts are a code
reference rather than a turnkey pipeline: input paths, caller output locations and several
reference files are written for the original project environment, so they need editing before the
analysis can be re-run elsewhere. The raw sequencing data are deposited separately, and the tables
each arm produces are included here, so the results can be inspected without re-running anything.

The figure package under `Figures/FigureFiles/` is the exception. It carries no absolute paths,
reads only files committed to this repository, and renders every manuscript figure from a fresh
clone with no configuration.

---

## Repository structure

```
LRS-FusionBench/
│
├── CellLineData/
│   Analysis of SGNex K562 and MCF7 cell line fusion datasets
│   SGNex_Cell_lines_K562_MCF7_analysis.rmd drives this arm
│
├── Figures/
│   ├── Input_dataframes/   every table the figure code reads
│   ├── FigureFiles/        the manuscript figure package (Python)
│   │   ├── code/           one folder per manuscript figure
│   │   └── output/         rendered SVG, named by manuscript figure number
│   └── SCRIPTS/            R scripts that write summary tables
│
├── Fusions/
│   Huh7 fusion preprocessing, discovery benchmarking, fusion labelling
│   and tool execution scripts
│
├── Known_Fusions/
│   DepMap fusion calls and curated Huh7 reference fusions
│
├── PreprocessingReads/
│   Shell scripts for adapter trimming, read filtering and preprocessing
│
├── ReadLengths/
│   Read length analyses and summary statistics
│
├── SimulatedData/
│   Scripts for analysing simulated benchmarking datasets
│   Data_Analysis_Scripts/simulated_reads_recall.rmd drives this arm
│
├── Preparing_data.rmd      drives the Huh7 arm
└── README.md
```

The analysis runs as three independent arms, each driven by one R notebook: the simulated
benchmark, the Huh7 libraries, and the SG-NEx cell lines. Each arm writes summary tables into
`Figures/Input_dataframes/`, which the figure code then reads.

---

## Repository contents

### PreprocessingReads/

Scripts for preparing sequencing reads prior to fusion detection.

* `DoradoTrimAdaptor.sh`: Adapter trimming using Dorado.
* `Fastp_trim_and_filter.sh`: Read trimming and quality filtering with Fastp.
* `Filter_NanoporeReads_Chopper.sh`: Read filtering using Chopper.
* `sequencing_data_shortcuts.sh`: Convenience commands for managing sequencing datasets.

---

### Fusions/

Contains the primary benchmarking workflows performed on our Huh7 Data.

Major components include:

* Formatting and standardisation of fusion caller outputs
* Fusion discovery benchmarking
* Fusion labelling and annotation
* Tool execution scripts

Key files:

* `FormattingandPreprocessFusions.R`
* `fusions_Huh7_discovery.tsv.gz`
* `fusions_readsupport_Huh7_discovery.tsv.gz`

Subdirectories include:

* `Discovery/`
* `label_fusions/`
* `runtools/`

---

### Known_Fusions/

Publicly available fusion calls together with manually curated Huh7 fusion references.

Included resources:

* `DepMap_OmicsFusionFiltered.csv` and `Model.csv`, from the DepMap 25Q2 public release
* `KnownHuh7FusionsLiterature_2026.csv`, Huh7 and HCC fusions curated from the literature
* `KnownHuh7Fusion_ah.csv`, the curated list annotated with Ensembl v110 identifiers
* `Cell_Line_Fusions/`, per cell line reference fusions

`ImportFormat_Known_Huh7.R` standardises the curated list and attaches Ensembl gene identifiers,
producing the reference set the benchmark matches against.

Redistribution terms for the DepMap release are recorded in `DATA_LICENSES.md`.

---

### CellLineData/

Analysis scripts for publicly available long-read sequencing datasets.

Includes analyses of:

* SG-NEx K562
* SG-NEx MCF7

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

**`Input_dataframes/`** holds every table the figure code reads. Most are written by the three
analysis arms; the remainder are read-level summaries computed from the alignments, a curated
performance table, and the reference annotation used to place genes on the ideogram.

**`SCRIPTS/`** holds the R scripts that compute summary tables during the simulated and SG-NEx
arms.

**`FigureFiles/`** is the manuscript figure package. `code/` has one folder per manuscript
figure, with a script per panel named after the panel it produces, so `code/Fig1/Fig1B.py`
renders Figure 1B. Shared modules sit alongside: `_paths.py` resolves every path relative to
itself, `palette.py` carries the shared style, and `_data.py`, `_known_fusion_data.py` and
`_similarity_helper.py` are shared loaders. `code/_make/` regenerates the two derived tables
that the figure code needs but the R pipeline does not write.

To render a figure, run its script from its own folder:

```bash
cd Figures/FigureFiles/code/Fig1 && python3 Fig1B.py
```

Output is written to `Figures/FigureFiles/output/<figure>/`. Paths are resolved relative to the
package, so the repository runs unchanged from any location. Vector SVG is tracked; the 600 dpi
PNG preview is regenerated on demand and is not.

Figure 1A and Supplementary Figure 14 are schematics drawn externally and have no code here.

---

## Figure source map

Each panel is rendered by the script named after it, under `Figures/FigureFiles/code/`. Panel
lettering is added when the panels are assembled, so it does not appear in the rendered files.

| Panel | Script | Input |
|---|---|---|
| Fig 1A | drawn externally | |
| Fig 1B | `Fig1/Fig1B.py` | `perf_value.tsv` |
| Fig 1C | `Fig1/Fig1C.py` | `stat_summary.tsv` |
| Fig 1D | `Fig1/Fig1D.py` | `stat_summary.tsv` |
| Fig 1E | `Fig1/Fig1E.py` | `fusion_profile_stat_summary.tsv` |
| Fig 1F | `Fig1/Fig1F.py` | `absolute_Breakpoint_Accuracy.tsv` |
| Fig 2A | `Fig2/Fig2A.py` | `fusions_in_n_algorithms.tsv` |
| Fig 2B | `Fig2/Fig2B.py` | `readsupp_filtering.tsv` |
| Fig 3A | `Fig3/Fig3A.py` | `counts_SGNex_summary.tsv.gz` |
| Fig 3B | `Fig3/Fig3B.py` | `fusionlevel_SGNex_summary.tsv.gz` |
| Fig 3C | `Fig3/Fig3C.py` | Huh7 discovery table |
| Fig 3D | `Fig3/Fig3D.py` | `call_burden.tsv` |
| Fig 3E | `Fig3/Fig3E_SFig12_circos.py` | Huh7 discovery table, `ref_annot.gtf.gene_spans` |
| Fig 3F | `Fig3/Fig3F.py` | `promiscuity_lib_metrics.tsv` |
| Fig 4A | `Fig4/Fig4A.py` | `huh7_readqc_hist_allreads.tsv` |
| Fig 4B | `Fig4/Fig4B.py` | `huh7_txcoverage_pertx.tsv` |
| Fig 4C | `Fig4/Fig4C.py` | Huh7 discovery table, `KnownHuh7Fusion_ah.csv` |
| Fig 4D | `Fig4/Fig4D.py` | Huh7 discovery table |
| Fig 4E | `Fig4/Fig4E.py` | Huh7 discovery table |
| Supp Fig 1 | `SFig1/SFig1.py` | `stat_summary.tsv` |
| Supp Fig 2 | `SFig2/SFig2.py` | `fusion_profile_stat_summary.tsv` |
| Supp Fig 3 | `SFig3/SFig3.py` | `stat_summary.tsv` |
| Supp Fig 4 | `SFig4/SFig4.py` | `fusion_profile_stat_summary.tsv` |
| Supp Fig 5A | `SFig5/SFig5A.py` | `fusion_profile_stat_summary.tsv` |
| Supp Fig 5B | `SFig5/SFig5B.py` | `fusion_profile_stat_summary.tsv` |
| Supp Fig 6 | `SFig6/SFig6.py` | `absolute_Breakpoint_Accuracy.tsv` |
| Supp Fig 7 | `SFig7/SFig7.py` | `fusions_in_n_algorithms.tsv` |
| Supp Fig 8A | `SFig8/SFig8A.py` | `readsupp_filtering.tsv` |
| Supp Fig 8B | `SFig8/SFig8B.py` | `readsupp_filtering.tsv` |
| Supp Fig 9 | `SFig9/SFig9.py` | Huh7 discovery table |
| Supp Fig 10 | `SFig10/SFig10.py` | Huh7 discovery table |
| Supp Fig 11A | `SFig11/SFig11A.py` | Huh7 discovery table |
| Supp Fig 11B | `SFig11/SFig11B.py` | Huh7 discovery table |
| Supp Fig 12 | `Fig3/Fig3E_SFig12_circos.py` | Huh7 discovery table, `ref_annot.gtf.gene_spans` |
| Supp Fig 13 | `SFig13/SFig13.py` | `huh7_read_identity_sample.tsv.gz` |
| Supp Fig 14 | drawn externally | |

One script renders both Figure 3E and Supplementary Figure 12: they are the same circos plot drawn
for every caller and library preparation, with the JAFFAL panels used in the main text and the
remaining five callers in the supplement.

Inputs named without a path are in `Figures/Input_dataframes/`; the Huh7 discovery table is
`Fusions/fusions_readsupport_Huh7_discovery.tsv.gz`.

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

The three analysis arms are written in **R**; the manuscript figures are written in **Python**.

Python packages used by the figure code:

* pandas
* numpy
* matplotlib
* scipy
* networkx
* pycirclize

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

The DepMap fusion release is redistributed for convenience. Its terms are recorded in `DATA_LICENSES.md`.

---

## Citation

If you use this repository in your work, please cite:

> Dorney R\*, Wu S\*, Hung JY-H, Hebbard L, Schmitz U. A comprehensive benchmark of
> transcriptome-wide fusion detection using long-read RNA sequencing. *bioRxiv* 2026.
> doi: [10.64898/2026.08.07.743439](https://doi.org/10.64898/2026.08.07.743439)

\* These authors contributed equally.

```bibtex
@article{Dorney2026LRFusionBenchmark,
  title   = {A comprehensive benchmark of transcriptome-wide fusion detection
             using long-read RNA sequencing},
  author  = {Dorney, Ryley and Wu, Siyuan and Hung, Julia Yun-Hsuan and
             Hebbard, Lionel and Schmitz, Ulf},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.08.07.743439}
}
```

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

