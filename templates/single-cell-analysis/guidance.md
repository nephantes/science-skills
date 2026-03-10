When analysing single-cell RNA-seq data:
- Load h5ad files with scanpy (sc.read_h5ad)
- Check adata.obs for cell type annotations and conditions
- Use sc.tl.rank_genes_groups for differential expression between conditions
- Filter DEGs by adjusted p-value < 0.05 and |log2FC| > threshold
- For pathway enrichment, use gseapy.enrichr or gseapy.prerank on DEG lists
- Report number of DEGs per cell type
- Use adata.raw if available for expression values (un-normalised counts)