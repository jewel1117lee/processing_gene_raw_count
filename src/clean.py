import pandas as pd
from openpyxl import Workbook

def clean_up(df, gene_name_col):
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

def filter_zeros(df: pd.DataFrame, subset_cols: list[str] | None = None) -> pd.DataFrame:
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
    

if __name__ == "__main__":
    gene_counts_path = "data/table_raw_count_for_all_samples_with_txID_gene_name.xlsx"
    protein_coding_path = "data/protein coding gene.tsv" 

    df1 = pd.read_excel(gene_counts_path)
    protein_coding = pd.read_csv(protein_coding_path, sep='\t')

    df1['gene_name'] = clean_up(df1['gene_name'])
    print(df1)

    df1 = df1[df1['gene_name'].isin(protein_coding['gene_name'])]
    print(df1)

    df1 = remove_all_zeros(df1, axis='columns')
    print(df1)

    output_path = "data/cleaned_counts.csv"
    df1.to_csv(output_path, index=False)