When performing correlation analysis, consider which method fits
the data distribution. Different methods capture different relationships.

METHOD SELECTION GUIDE:
- Pearson: For linear relationships between normally distributed variables. Sensitive to outliers.
- Spearman: For monotonic relationships. Robust to outliers and non-normality. Use when data has outliers or is ordinal.
- Kendall: For small samples or many ties. More robust than Spearman but slower.

ALWAYS:
- Check for and handle outliers before Pearson (or use Spearman instead)
- For large feature matrices (>5000 features), subset to top-variance features to avoid memory issues
- Report both correlation coefficient AND p-value
- Use clustermap (hierarchical-clustered heatmap) for visualization
- Save the correlation matrix to CSV