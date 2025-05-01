# Processing Raw Counts

A simple end-to-end pipeline for RNA-seq count processing, normalization, differential expression, and downstream comparison to published gene lists. Built in Python, it organizes reusable modules under `src/` and provides a single `main.py` driver.

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
  - main.py ← entry-point driver script

## Usage
I ran it with python main.py in terminal, but original code is in jupyter notebook which is a lot easier to use for bioinformatics with a lot of calling the functions with different data and graphing

## Pipeline Overview
The main.py has a much better explanation, but overall my process is 
1. Read in the data
2. Clean the data
    2a. filter noncoding genes
3. Produce RPKM
    3a. clean the data again
    3b. run the function for RPKM
    3c. export
4. Produce DEG analysis
    4a. clean the data/ normalization
    4b. run the function for DEG analysis
    4c. export 
5. Compare to data publish online
    5a. clean the data that is published online
    5b. fun the function to compare
    5c. run statistical test, such as t-test
    5d. graph
    5e. export

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



