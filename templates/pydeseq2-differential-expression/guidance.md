When performing differential expression in Python, use pydeseq2 (DESeq2 algorithm).
Always do all pairwise comparisons via itertools.combinations(). When there are >1 comparison,
ALWAYS create the DEG heatmap (deg_heatmap.png) and direction barplot (deg_barplot.png) BEFORE
any volcano plots. Sample metadata is loaded from sample_metadata.csv; the first metadata column
is used as the group variable. Round counts to integers before passing to DeseqDataSet.
