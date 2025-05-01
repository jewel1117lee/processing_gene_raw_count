# analysis.py: Pipeline to match gene data with published data, perform t-test on matched vs unmatched genes,
# plot results, and export zipped data.

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
from openpyxl import Workbook
from clean import clean_up, filter_noncoding, filter_zeros, dropna_from_lists


def match_database(df, gene_names, values, comp_df, comp_col, protein_coding_df):
    """
    Match genes in gene_names to comp_df[comp_col] after filtering protein-coding genes.

    Args:
        gene_names (list of str): List of gene names.
        values (list of numeric): Corresponding values for each gene (e.g., log-RPKM).
        comp_df (pd.DataFrame): Published data containing a column comp_col with gene names.
        comp_col (str): Column name in comp_df to match against gene_names.
        protein_coding_df (pd.DataFrame): DataFrame with 'gene_name' column listing protein-coding genes.

    Returns:
        tuple: (matched_genes, matched_values, unmatched_genes, unmatched_values, non_matched_in_pub)
    """
    # Prepare published gene list (uppercase)
    pub_genes = clean_up(comp_df,comp_col)[comp_col]
    
    df_new = clean_up(df, gene_names)

    df_gene_value = df_new[[gene_names, values]]

    mask = df_gene_value[gene_names].isin(pub_genes)
    
    matched   = df_gene_value[mask]
    unmatched = df_gene_value[~mask]

    matched_genes = matched[gene_names].tolist()
    matched_vals = matched[values].tolist()
    unmatched_genes = unmatched[gene_names].tolist()
    unmatched_vals = unmatched[values].tolist()

    return matched_genes, matched_vals, unmatched_genes, unmatched_vals


def ttest(data1, data2):
    """
    Perform Welch's t-test between two groups and adjust p-value by factor of number of tests.

    Args:
        data1, data2 (list or array-like): Numeric values for two groups.

    Returns:
        tuple: ((t_stat, p_adj), significance_star)
    """
    t_stat, p_val = scipy.stats.ttest_ind(data1, data2, equal_var=False, nan_policy="omit")
    # Example p-value adjustment (Bonferroni by 2 tests)
    p_adj = p_val 

    print()
    if np.isnan(data1).any() or np.isnan(data2).any():
        print("There is NaN value in the data you sent")
    
    # Significance stars
    if p_adj < 0.0001:
        star = '****'
    elif p_adj < 0.001:
        star = '***'
    elif p_adj < 0.01:
        star = '**'
    elif p_adj < 0.05:
        star = '*'
    else:
        star = ''

    return (t_stat, p_adj, star)


def make_box_figure(
    title, data, labels, xlabel, ylabel,
    colors=None, median_color='black', ttest= None, y = None, output_path=None
):
    """
    Create and optionally save a boxplot comparing two datasets.

    Args:
        title (str): Plot title.
        data (list of lists): Two lists of numeric values.
        labels (list of str): Labels for the two boxes.
        xlabel, ylabel (str): Axis labels.
        colors (list of str): Face colors for the boxes.
        median_color (str): Color for median lines.
        ttest: Any statistical test that you will like above your box blot.
        output_path (str): If provided, saves plot to this path.

    Returns:
        matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    data = dropna_from_lists(data)
    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(
        data,
        showfliers=False,
        patch_artist=True
    )
    # Apply colors if provided
    if colors:
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
    # Median line color
    for median in bp['medians']:
        median.set_color(median_color)

    ax.set(title = title, xlabel= xlabel, ylabel = ylabel)
    ax.set_xticklabels(labels)
    ax.yaxis.grid(True, linestyle='-', color='lightgrey', alpha=0.5)

    if y is None:
        # compute the 75th percentile for each group
        q3s = [np.percentile(group, 75) for group in data]
        # take the highest one
        max_q3 = max(q3s)
        # then position your annotation relative to that
        y = max_q3 * 4.5
            
    for i, ttest in enumerate(ttest, start=1):
        if type(ttest) == str: 
            ax.text(i, y, ttest, ha='center', va='bottom', color='black')
        if type(ttest) == float: 
            ax.text(i, y, f"{ttest:.2e}", ha='center', va='bottom', color='black')
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.savefig(output_path)

    return fig, ax


def zip_data(data, col_name=None, output_file=None):  
    ""
    max_length = max(len(lst) for lst in data)
    for lst in data:
        while len(lst) < max_length:
            lst.append(None)
    
    data_all = pd.DataFrame(data)
    data_all = data_all.transpose()

    if col_name:
        if len(col_name) != len(data_all.columns):
            raise ValueError("name column must be equal to data columns")
        else:
            data_all.columns = col_name

    if output_file:
        data_all.to_excel(output_file, index=False)

    return(data_all)

if __name__ == "__main__":
    main()