# Processing Raw Counts

A simple end-to-end pipeline for RNA-seq count processing, normalization, differential expression, and downstream comparison to published gene lists. Built in Python, it organizes reusable modules under `src/` and provides a single `main.py` driver.

---

## Table of Contents

1. [Features](#features)  
2. [Directory Structure](#directory-structure)  
3. [Usage](#usage)  
4. [Pipeline Overview](#pipeline-overview)  
5. [Configuration & Paths](#configuration--paths)  
6. [Module Reference](#module-reference)  
7. [Contributing](#contributing)  
8. [License](#license)

---

## Features

- **Filter non-coding & zero counts**  
- **Compute RPKM & CPM** normalization  
- **Run differential expression** with PyDESeq2  
- **Match results** against external gene lists  
- **Statistical testing & visualization** (t-test + boxplots)  

---

## Directory Structure

processing_raw_count/ ├── data/
│ ├── table_raw_count_for_all_samples_with_txID_gene_name.xlsx
│ ├── protein coding gene.tsv
│ └── mouse_table_DF.csv
├── result/ ← output folder for Excel & plots
├── src/
│ ├── init.py
│ ├── clean.py ← filter_noncoding, filter_zeros, dropna_from_lists
│ ├── normalization.py ← compute_rpkm, compute_cpm
│ ├── analysis.py ← match_database, ttest, make_box_figure, zip_data
│ └── deg_analysis.py ← run_deg_analysis
└── main.py ← entry-point driver script
