from typing import List, Optional, Union
import numpy as np
import pandas as pd

def calc_rpkm(read_count, gene_len, total_reads):
    """
    Calculate RPKM (Reads Per Kilobase per Million mapped reads).

    Args:
        read_count: Number of reads mapped to a gene (scalar or array-like).
        gene_len:   Length of the gene in base-pairs (scalar or array-like).
        total_reads: Total number of mapped reads in the experiment (scalar).

    Returns:
        RPKM value(s), same type as inputs (float, numpy array, or pandas Series).
    """
    # Prevent division by zero
    if np.any(gene_len == 0):
        raise ValueError("gene_len contains zero(s), cannot divide by zero.")
    if total_reads == 0:
        raise ValueError("total_reads must be > 0.")

    # Convert everything to float (handles both scalars and arrays)
    gene_kb = gene_len / 1_000.0
    million_reads = total_reads / 1_000_000.0

    return read_count / (gene_kb * million_reads)


def compute_rpkm(
    df: pd.DataFrame,
    gene_length_col: str = 'gene_length',
    sample_cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Compute RPKM for each sample column in `df`.

    Args:
        df: DataFrame containing raw counts and a gene-length column.
        gene_length_col: Name of the column with gene lengths (bp).
        sample_cols: List of columns to treat as count data. 
                     If None, all numeric cols except `gene_length_col` are used.

    Returns:
        A new DataFrame with additional columns "<sample>_RPKM" for each sample.
    """
    # Copy to avoid mutating original
    df_out = df.copy()
    
    # Ensure gene_length exists
    if gene_length_col not in df_out:
        raise KeyError(f"Column '{gene_length_col}' not found in DataFrame.")
    
    # Determine which columns to process
    if sample_cols is None:
        # pick numeric columns except gene_length_col
        sample_cols = [
            c for c in df_out.columns
            if c != gene_length_col and pd.api.types.is_numeric_dtype(df_out[c])
        ]
    else:
        # validate user-specified columns
        missing = [c for c in sample_cols if c not in df_out.columns]
        if missing:
            raise KeyError(f"Sample columns not found: {missing}")
    
    # vector of gene lengths
    gene_len = df_out[gene_length_col].astype(float)
    
    # compute RPKM per sample
    for col in sample_cols:
        total = df_out[col].sum()
        df_out[f"{col}_RPKM"] = calc_rpkm(df_out[col].astype(float), gene_len, total)
    
    return df_out


def calc_cpm(
    read_count: Union[pd.Series, np.ndarray, float],
    total_reads: float):
    """
    Calculate CPM (Counts Per Million) for a sample or vector of counts.

    Args:
        read_count: Number of reads mapped to a feature (scalar or array-like).
        total_reads: Total number of mapped reads in the sample (scalar).

    Returns:
        CPM value(s), same type as read_count (float, np.ndarray, or pd.Series).

    Raises:
        ValueError: if total_reads <= 0.
    """
    if total_reads <= 0:
        raise ValueError("total_reads must be > 0")

    factor = 1_000_000 / total_reads
    return read_count * factor

def compute_cpm_df(
    df: pd.DataFrame,
    sample_cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Compute CPM for each sample column in `df`.

    Args:
        df: DataFrame containing raw count columns.
        sample_cols: List of columns to treat as count data.
                     If None, auto-detects all numeric columns.

    Returns:
        A new DataFrame with one "<sample>_CPM" column per sample.
    """
    # Determine which columns are counts
    if sample_cols is None:
        sample_cols = [
            c for c in df.columns
            if pd.api.types.is_numeric_dtype(df[c])
        ]
    else:
        missing = [c for c in sample_cols if c not in df.columns]
        if missing:
            raise KeyError(f"Sample columns not found: {missing}")

    df_out = pd.DataFrame(index=df.index)
    for col in sample_cols:
        total = float(df[col].sum())
        df_out[f"{col}_CPM"] = calc_cpm(df[col].astype(float), total)

    return df_out
