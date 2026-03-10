import pandas as pd
from scipy.stats import spearmanr
import gseapy as gp

# Load MAGeCK results
data = pd.read_excel("mageck_results.xlsx", sheet_name=0)
print(f"Total genes: {len(data)}")

# Rank genes by score
ranked = data.sort_values("neg|score", ascending=False)
sig_genes = ranked[ranked["neg|fdr"] < 0.05]
print(f"Significant hits (FDR < 0.05): {len(sig_genes)}")

# Spearman correlation between conditions
rho, pval = spearmanr(data["condition1_score"], data["condition2_score"])
print(f"Spearman rho: {rho:.4f}, p-value: {pval:.2e}")

# Pathway enrichment
gene_list = sig_genes["gene"].tolist()
enr = gp.enrichr(gene_list=gene_list, gene_sets="Reactome_2022", organism="human")
enr.results.to_csv("pathway_enrichment.csv", index=False)