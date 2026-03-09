When performing clustering analysis, use multiple methods and compare
results. Clusters that appear consistently across methods are more reliable.

METHOD SELECTION GUIDE:
- Hierarchical clustering: Best for visualizing relationships (dendrograms). Use ward linkage as default.
- K-means: Best for large datasets when you have a target k. Always run silhouette analysis.
- DBSCAN: Best when clusters have irregular shapes or you don't know k. Requires epsilon tuning.

ALWAYS:
- Standardize/scale data before clustering (StandardScaler)
- Remove zero-variance features first
- Try at least 2 methods and compare
- Use silhouette score or similar metric to evaluate cluster quality
- Visualize clusters in 2D (PCA or UMAP space)
- Save cluster assignments to CSV for downstream use