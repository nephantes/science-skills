IMPORTANT: First check if the uploaded file already contains pre-computed DE results
(columns like log2FoldChange/logFC and padj/FDR/adj.P.Val). If so, skip pydeseq2 entirely
and proceed directly to downstream analysis using those columns.

When performing differential expression in Python, use pydeseq2. Ensure sample metadata
is loaded from sample_metadata.csv and contrasts are derived from the data, never hardcoded.
The count matrix must have genes as rows and samples as columns. Round counts to integers
before passing to DeseqDataSet.

MULTIPLE COMPARISONS IN INPUT: If the input already has multi-comparison results (a 'comparison'
column or wide-format prefixed columns), loop over each comparison to generate per-comparison
volcano/MA plots, naming each file with the comparison label.