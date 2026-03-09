import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

expr = pd.read_csv("expression_matrix.csv", index_col=0)
print(f"Input: {expr.shape[0]} features x {expr.shape[1]} samples", flush=True)

# Limit to top-variance features for performance
if expr.shape[0] > 2000:
    var = expr.var(axis=1)
    expr = expr.loc[var.nlargest(2000).index]
    print(f"Subsetted to top 2000 variable features", flush=True)

# Method 1: Pearson correlation (sample-sample)
corr_pearson = expr.corr(method='pearson')
corr_pearson.to_csv('correlation_pearson.csv')
print(f"Pearson range: [{corr_pearson.min().min():.3f}, {corr_pearson.max().max():.3f}]", flush=True)

# Method 2: Spearman correlation (robust to outliers)
corr_spearman = expr.corr(method='spearman')
corr_spearman.to_csv('correlation_spearman.csv')
print(f"Spearman range: [{corr_spearman.min().min():.3f}, {corr_spearman.max().max():.3f}]", flush=True)

# Compare methods
diff = (corr_pearson - corr_spearman).abs()
print(f"Mean abs difference Pearson vs Spearman: {diff.mean().mean():.4f}", flush=True)
print(f"Max abs difference: {diff.max().max():.4f}", flush=True)

# Visualize with clustermap (Spearman — more robust)
g = sns.clustermap(corr_spearman, cmap='RdBu_r', vmin=-1, vmax=1,
                   figsize=(10, 10), linewidths=0)
plt.savefig("correlation_clustermap.png", dpi=100, bbox_inches="tight")
plt.close()