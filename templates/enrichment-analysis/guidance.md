When performing functional enrichment analysis, choose the method
best suited to the input. Present results from at least 2 approaches when possible.

METHOD SELECTION GUIDE:
- Over-representation analysis (ORA): Best when you have a discrete gene list (e.g., DE genes with padj < 0.05). Use scipy.stats.fisher_exact or gseapy.enrichr.
- Gene Set Enrichment Analysis (GSEA): Best when you have ranked data (e.g., all genes ranked by fold change). Use gseapy.prerank.
- If the gene list is very small (<20 genes), prefer ORA with Fisher's exact test over GSEA.
- If the gene list is large (>500 genes), consider tightening the significance threshold first.

ALWAYS:
- Report both enrichment score/odds ratio AND adjusted p-value
- Use multiple testing correction (Benjamini-Hochberg)
- Show top 15-20 enriched terms in a bar plot or dot plot
- Save full results table to CSV