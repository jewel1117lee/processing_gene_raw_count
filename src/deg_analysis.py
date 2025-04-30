import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from typing import List, Optional, Union


def run_deg_analysis(
    counts_df: pd.DataFrame,
    group1_cols: List[str],
    group2_cols: List[str],
    g1_label: str = "Group1",
    g2_label: str = "Group2",
    gene_names: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Differential expression (GLM-based) comparing two sets of samples.

    Args:
        counts_df: DataFrame with gene counts. Rows=genes, cols=samples.
        group1_cols: List of column names in counts_df for group 1.
        group2_cols: List of column names in counts_df for group 2.
        g1_label: Name to assign to group1 (for design matrix).
        g2_label: Name to assign to group2.
        gene_names: Optional list of gene names (len==n_genes). 
                    If None, uses counts_df.index as gene names.

    Returns:
        DataFrame with columns ['gene', 'logFC', 'pval', 'FDR'], sorted by FDR.
    """
    # === 1. Subset and validate ===
    all_cols = list(group1_cols) + list(group2_cols)
    missing = [c for c in all_cols if c not in counts_df.columns]
    if missing:
        raise KeyError(f"These samples are missing in counts_df: {missing}")

    sub = counts_df[all_cols]
    n1, n2 = len(group1_cols), len(group2_cols)

    # === 2. Prepare gene names ===
    if gene_names is not None:
        if len(gene_names) != sub.shape[0]:
            raise ValueError("Length of gene_names must match number of rows in counts_df.")
        genes = gene_names
    else:
        genes = list(sub.index)

    # === 3. Build AnnData for normalization + log1p ===
    adata = sc.AnnData(X=sub.transpose())
    adata.var_names = genes
    adata.obs['group'] = [g1_label]*n1 + [g2_label]*n2

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # === 4. Design matrix: intercept + one dummy for group2 vs group1 ===
    design = pd.get_dummies(adata.obs['group'], drop_first=True)
    design = sm.add_constant(design)

    # === 5. Fit GLM per gene ===
    results = []
    for gene in adata.var_names:
        # y = normalized counts (log1p) per sample
        mat = adata[:, gene].X
        y = mat.toarray().flatten() if sparse.issparse(mat) else mat.flatten()

        model = sm.GLM(y, design, family=sm.families.Poisson())
        fit = model.fit()

        # coef[1] is the effect of the second dummy column (group2 vs group1)
        coef = fit.params[1]
        pval = fit.pvalues[1]
        results.append((gene, coef, pval))

    # === 6. Assemble DataFrame, adjust p-values ===
    res_df = pd.DataFrame(results, columns=["gene", "logFC", "pval"])
    res_df["FDR"] = multipletests(res_df["pval"], method="fdr_bh")[1]
    res_df = res_df.sort_values("FDR").reset_index(drop=True)

    return res_df