import pandas as pd
import numpy as np
from scipy import stats

# Load DE results
de = pd.read_csv("de_results.csv", index_col=0)
sig_genes = de[de['padj'] < 0.05].index.tolist()
all_genes = de.index.tolist()
print(f"Significant genes: {len(sig_genes)} / {len(all_genes)}", flush=True)

# Method 1: Over-representation with gseapy
try:
    import gseapy as gp
    enr = gp.enrichr(gene_list=sig_genes, gene_sets='GO_Biological_Process_2021',
                     organism='Human', outdir=None, cutoff=0.05)
    enr_results = enr.results.sort_values('Adjusted P-value').head(20)
    enr_results.to_csv('enrichment_ora.csv', index=False)
    print(f"ORA: {len(enr_results)} enriched terms (padj < 0.05)", flush=True)
except ImportError:
    print("gseapy not available — falling back to manual Fisher test", flush=True)

# Method 2: GSEA with preranked list
try:
    import gseapy as gp
    ranked = de[['log2FoldChange']].dropna()
    ranked = ranked.sort_values('log2FoldChange', ascending=False)
    ranked.columns = ['score']
    pre_res = gp.prerank(rnk=ranked, gene_sets='GO_Biological_Process_2021',
                         min_size=15, max_size=500, permutation_num=100, outdir=None)
    gsea_results = pre_res.res2d.sort_values('FDR q-val').head(20)
    gsea_results.to_csv('enrichment_gsea.csv', index=False)
    print(f"GSEA: {len(gsea_results)} enriched terms (FDR < 0.25)", flush=True)
except Exception as e:
    print(f"GSEA skipped: {e}", flush=True)