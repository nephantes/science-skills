When performing differential expression in Python, use pydeseq2.
Ensure sample metadata is loaded from sample_metadata.csv and contrasts are
derived from the data, never hardcoded. The count matrix must have genes as
rows and samples as columns. Round counts to integers before passing to DeseqDataSet.