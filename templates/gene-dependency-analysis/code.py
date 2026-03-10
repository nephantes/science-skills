import pandas as pd
import numpy as np
from scipy.stats import spearmanr, skew
from statsmodels.stats.multitest import multipletests

# Load data
expression = pd.read_csv("expression.csv", index_col=0)
dependency = pd.read_csv("gene_effect.csv", index_col=0)

# Align samples
common_lines = list(set(expression.index) & set(dependency.index))
expr = expression.loc[common_lines]
dep = dependency.loc[common_lines]

# Distribution characterisation
expr_skew = skew(expr.values.flatten(), nan_policy="omit")
print(f"Expression skewness: {expr_skew:.2f} (right-skewed)" if expr_skew > 0 else f"Expression skewness: {expr_skew:.2f}")

# Compute per-gene correlations
results = []
common_genes = list(set(expr.columns) & set(dep.columns))
for gene in common_genes:
    mask = expr[gene].notna() & dep[gene].notna()
    if mask.sum() < 10:
        continue
    rho, pval = spearmanr(expr.loc[mask, gene], dep.loc[mask, gene])
    results.append({"gene": gene, "spearman_rho": rho, "pvalue": pval})

res_df = pd.DataFrame(results)
_, res_df["fdr"], _, _ = multipletests(res_df["pvalue"], method="fdr_bh")
res_df = res_df.sort_values("spearman_rho")
res_df.to_csv("gene_dependency_correlations.csv", index=False)
print(f"Strongest negative: {res_df.iloc[0]['gene']} (rho={res_df.iloc[0]['spearman_rho']:.4f})")