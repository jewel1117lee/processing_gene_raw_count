import pandas as pd
from typing import List
import pydeseq2.preprocessing as preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_deg_analysis(
    counts_df: pd.DataFrame,
    gene_id_col:str,
    group1_cols: List[str],
    group2_cols: List[str],
    g1_label: str = "Group1",
    g2_label: str = "Group2",
    gene_name_col: str = 'gene_name',
    normalize_factor = False
) -> pd.DataFrame:
    """
    Differential expression between two sets of columns using PyDESeq2.

    Args:
        counts_df:    DataFrame of raw counts (genes × samples).
        group1_cols:  List of sample‐column names for group 1.
        group2_cols:  List of sample‐column names for group 2.
        g1_label:     Label to assign to group1 in the design.
        g2_label:     Label to assign to group2.

    Returns:
        A DataFrame of DESeq2‐style results (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj).
    """
    gene_name_df = counts_df[[gene_id_col,gene_name_col]]
    counts_df = counts_df.set_index(gene_id_col)
    
    # 1) Check columns exist
    all_cols = group1_cols + group2_cols

    # 2) Subset to only those columns
    sub_counts = counts_df.loc[:, all_cols]
    counts = sub_counts.T

    # 3) Build a small metadata (clinical) table
    metadata = pd.DataFrame({
        "condition": [g1_label] * len(group1_cols) + [g2_label] * len(group2_cols)
    }, index=all_cols)

    if normalize_factor == True: 
        norm_counts, lib_sizes = preprocessing.deseq2_norm(counts)
        norm_counts = norm_counts.round().astype(int)
        dds = DeseqDataSet(counts=norm_counts, metadata=metadata, design_factors="condition")
    
    else:
        dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition")
        
    dds.deseq2()
    dds.fit_LFC()
    
    stat_res = DeseqStats(dds, contrast = ('condition', g1_label, g2_label))
    stat_res.summary()
    
    res = stat_res.results_df

    res_df = gene_name_df.merge(res, on= gene_id_col, how="left")
    return(res_df)

if __name__ == "__main__":
    main()