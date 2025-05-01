# Processing Raw Counts

A simple end-to-end pipeline for RNA-seq count processing, normalization, differential expression, and downstream comparison to published gene lists. Built in Python, it organizes reusable modules under `src/` and provides a single `main.py` as a EXAMPLE driver.

---

## Table of Contents

1. [Features](#features)  
2. [Directory Structure](#directory-structure)  
3. [Usage](#usage)  
4. [Pipeline Overview](#pipeline-overview)  
5. [Module Reference](#module-reference)  

---

## Features

- **Filter non-coding & zero counts**  
- **Compute RPKM & CPM** normalization  
- **Run differential expression** with PyDESeq2  
- **Match results** against external gene lists  
- **Statistical testing & visualization** (t-test + boxplots)  

---

## Directory Structure

- processing_raw_count/
  - data/
    - raw_count_file
    - protein_coding_file (optional)
    - transcript_length_file (optional)
  - result/ ← output folder for Excel & plots
  - src/
    - __init__.py
    - clean.py ← filter_noncoding, filter_zeros, dropna_from_lists
    - normalization.py ← compute_rpkm, compute_cpm
    - analysis.py ← match_database, ttest, make_box_figure, zip_data
    - deg_analysis.py ← run_deg_analysis
  - main.py ← entry-point driver script (*** AN EXAMPLE ***)

## Usage
I ran it with python main.py in terminal, but original code is in jupyter notebook which is a lot easier to use for bioinformatics with a lot of calling the functions with different data and graphing

## Pipeline Overview
1. Read in the data  
2. Clean the data  
   1. Filter noncoding genes  
3. Produce RPKM  
   1. Clean the data again  
   2. Run the RPKM function  
   3. Export results  
4. Produce DEG analysis  
   1. Clean and normalize the data  
   2. Run the DEG analysis function  
   3. Export results  
5. Compare to published data  
   1. Clean the published data  
   2. Run the comparison function  
   3. Run statistical test (e.g., t-test)  
   4. Generate graphs  
   5. Export comparison results  

## Module Reference
src/clean.py

filter_noncoding(df, prot_df, gene_col, external_col)

filter_zeros(df)

dropna_from_lists(list_of_lists)

src/normalization.py

compute_rpkm(df, gene_id_cols, gene_length_col, sample_cols=None)

compute_cpm(df, sample_cols=None)

src/deg_analysis.py

run_deg_analysis(counts_df, gene_id_col, group1_cols, group2_cols, g1_label, g2_label, normalize_factor=False)

src/analysis.py

match_database(df, gene_names, values, comp_df, comp_col, protein_coding_df)

ttest(array1, array2)

make_box_figure(...)

zip_data(data, col_names=None)



