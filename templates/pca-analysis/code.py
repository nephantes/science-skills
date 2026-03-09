import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

expr = pd.read_csv("expression_matrix.csv", index_col=0)
meta = pd.read_csv("sample_metadata.csv")

X = expr.values.T  # samples as rows
X = X[:, X.var(axis=0) > 0]  # remove zero-variance
assert np.isfinite(X).all(), "NaN/inf in input"

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

n_comp = min(2, X.shape[0] - 1, X.shape[1])
pca = PCA(n_components=n_comp)
pcs = pca.fit_transform(X_scaled)

fig, ax = plt.subplots(figsize=(8, 6))
for grp in meta['condition'].unique():
    mask = meta['condition'].values == grp
    ax.scatter(pcs[mask, 0], pcs[mask, 1], label=grp, s=60, alpha=0.8)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.legend()
plt.savefig("pca_plot.png", dpi=100, bbox_inches="tight")
plt.close()