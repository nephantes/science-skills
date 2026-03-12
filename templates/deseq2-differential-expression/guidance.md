IMPORTANT: First check if the uploaded file already contains pre-computed DE results
(columns like log2FoldChange/logFC and padj/FDR). If so, skip DESeq2 entirely and proceed
directly to downstream analysis using those columns.

When performing differential expression with DESeq2, follow this tested workflow.
Ensure sample metadata is loaded from sample_metadata.csv and contrasts are derived from
the data, never hardcoded.

MULTIPLE COMPARISONS IN INPUT: If the input already has multi-comparison results (a 'comparison'
column or wide-format prefixed columns), loop over each comparison to generate per-comparison
volcano/MA plots, naming each file with the comparison label.