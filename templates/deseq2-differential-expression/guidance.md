When performing differential expression with DESeq2 in R, use this multi-group workflow.
Always do all pairwise comparisons via combn(). When there are >1 comparison, ALWAYS create
the DEG heatmap (deg_heatmap.png) and direction barplot (deg_barplot.png) BEFORE any volcano plots.
Sample metadata is loaded from sample_metadata.csv; the first non-'sample' column is used as the group variable.