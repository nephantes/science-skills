import scanpy as sc
import pandas as pd
import numpy as np

# Load the h5ad file
adata = sc.read_h5ad("data.h5ad")
print(f"Shape: {adata.shape}")
print(f"Cell types: {adata.obs.columns.tolist()}")

# Differential expression per cell type
sc.tl.rank_genes_groups(adata, groupby="condition", method="wilcoxon")
result = sc.get.rank_genes_groups_df(adata, group=None)
sig = result[result["pvals_adj"] < 0.05]
print(f"Significant DEGs: {len(sig)}")
sig.to_csv("de_results.csv", index=False)

# Optional: pathway enrichment with gseapy
import gseapy as gp
gene_list = sig["names"].tolist()
enr = gp.enrichr(gene_list=gene_list, gene_sets="GO_Biological_Process_2021", organism="human")
enr.results.to_csv("enrichment_results.csv", index=False)