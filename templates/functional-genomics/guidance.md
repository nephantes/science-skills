When analysing CRISPR screen / functional genomics data:
- MAGeCK output files contain gene-level scores (pos|neg|fdr columns)
- Load MAGeCK results from Excel/CSV and sort by score or FDR
- For correlation analysis: use scipy.stats.spearmanr or pearsonr
- For pathway enrichment: use gseapy prerank or enrichr with Reactome/KEGG gene sets
- GMT files can be loaded directly for custom gene set collections
- Use mygene for gene ID conversion if needed (pip install mygene)
- Report overlap between conditions, number of significant hits