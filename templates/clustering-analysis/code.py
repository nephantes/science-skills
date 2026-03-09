import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

expr = pd.read_csv("expression_matrix.csv", index_col=0)
X = expr.values.T  # samples as rows
X = X[:, X.var(axis=0) > 0]  # remove zero-variance features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
assert np.isfinite(X_scaled).all(), "NaN/inf after scaling"
print(f"Clustering matrix: {X_scaled.shape[0]} samples x {X_scaled.shape[1]} features", flush=True)

# Method 1: Hierarchical clustering
Z = linkage(X_scaled, method='ward')
fig, ax = plt.subplots(figsize=(12, 5))
dendrogram(Z, labels=expr.columns.tolist(), leaf_rotation=90, ax=ax)
plt.title("Hierarchical Clustering (Ward)")
plt.tight_layout()
plt.savefig("dendrogram.png", dpi=100, bbox_inches="tight")
plt.close()

# Method 2: K-means with silhouette analysis
sil_scores = []
K_range = range(2, min(10, X_scaled.shape[0]))
for k in K_range:
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = km.fit_predict(X_scaled)
    sil_scores.append(silhouette_score(X_scaled, labels))
best_k = list(K_range)[np.argmax(sil_scores)]
print(f"Best k={best_k} (silhouette={max(sil_scores):.3f})", flush=True)

km_final = KMeans(n_clusters=best_k, random_state=42, n_init=10)
clusters = km_final.fit_predict(X_scaled)
assignments = pd.DataFrame({'sample': expr.columns, 'cluster': clusters})
assignments.to_csv('cluster_assignments.csv', index=False)
print(f"Cluster sizes: {assignments['cluster'].value_counts().to_dict()}", flush=True)