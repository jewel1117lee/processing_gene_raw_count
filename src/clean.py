import pandas as pd
from typing import List, Optional
from openpyxl import Workbook

def capitalize_gene(df, gene_name_col):
    """
    """
    df = df.copy()
    df[gene_name_col] = df[gene_name_col].apply(lambda x: x.upper() if isinstance(x, str) else x)
    return df

def filter_noncoding(df, protein_coding, gene_name_col, protein_gene_name_col):
    """
    Return only the rows of df whose gene_name (after uppercasing)
    is found in the protein_coding list/DF. Everything else,
    including NaN or non-strings, is filtered out.
    """
    # Build the uppercase set of protein-coding names
    if isinstance(protein_coding, pd.DataFrame):
        prot_set = set(protein_coding[protein_gene_name_col].str.upper())
    else:
        prot_set = {str(g).upper() for g in protein_coding}

    # Define mask: True only if x is a string and its uppercase is in prot_set
    mask = df[gene_name_col].apply(
        lambda x: isinstance(x, str) and x.upper() in prot_set
    )

    # Return only matching rows
    return df[mask].copy()

def filter_zeros(df: pd.DataFrame, subset_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Remove rows from df that contain a zero in any of the specified columns.

    Args:
        df: Input DataFrame.
        subset_cols: List of column names to check for zeros. 
                     If None, all columns are checked.

    Returns:
        A new DataFrame with rows containing zeros in the checked columns removed.
    """
    # Work on a copy
    df_clean = df.copy()

    # Determine which columns to inspect
    if subset_cols is None:
        cols_to_check = df_clean.columns
        
    else:
        # Only keep those that actually exist in df
        cols_to_check = [c for c in subset_cols if c in df_clean.columns]
        if not cols_to_check:
            raise ValueError("No valid columns to check for zeros.")

    # Build mask: True for rows where *none* of the checked cols are zero
    mask = df_clean[cols_to_check].ne(0).all(axis=1)

    return df_clean.loc[mask].reset_index(drop=True)
    
def dropna_from_lists(data):
    """
    Given a list of lists, return a new list of lists
    where all `NaN` values have been removed.
    """
    return [
        [x for x in subgroup if x == x]   # `NaN != NaN`, so this drops NaNs
        for subgroup in data
    ]
