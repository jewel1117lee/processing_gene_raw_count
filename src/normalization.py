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
    gene_id_cols: Union[str, List[str]],
    gene_length_col: str = 'gene_length',
    sample_cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Returns a new DataFrame with:
      - the gene identifier column(s)
      - one <sample>_RPKM column per sample in sample_cols (or auto-detected)

    Args:
        df: Original DataFrame with counts and lengths.
        gene_id_cols: Column name or list of column names to carry over (e.g. 'gene_id').
        gene_length_col: Column holding gene lengths (bp).
        sample_cols: List of count columns. If None, auto-uses all numeric columns
                     except gene_length_col and gene_id_cols.

    Returns:
        New DataFrame with only the ID column(s) and the RPKM columns.
    """
    # Normalize inputs
    if isinstance(gene_id_cols, str):
        gene_id_cols = [gene_id_cols]
    
    # Detect sample columns if not provided
    if sample_cols is None:
        sample_cols = [
            c for c in df.columns 
            if c not in gene_id_cols + [gene_length_col]
            and pd.api.types.is_numeric_dtype(df[c])
        ]
    
    # Prepare output
    out = df[gene_id_cols].copy()
    gene_len = df[gene_length_col].astype(float)
    
    # Compute RPKM per sample and add to out
    for col in sample_cols:
        total = df[col].sum()
        rpkm = calc_rpkm(df[col].astype(float), gene_len, total)
        out[f"{col}_RPKM"] = rpkm
    
    return out


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

def compute_cpm(
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

if __name__ == "__main__":
    main()