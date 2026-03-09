library(tidyverse)
library(clusterProfiler)

# Load DESeq2 results (RDS or CSV)
de_results <- readRDS("de_results.rds")
de_df <- as.data.frame(de_results)

# Extract significant genes
sig_up <- de_df %>%
  filter(log2FoldChange > 1.5, padj < 0.05) %>%
  rownames_to_column("GeneID")
sig_dn <- de_df %>%
  filter(log2FoldChange < -1.5, padj < 0.05) %>%
  rownames_to_column("GeneID")

cat("Upregulated:", nrow(sig_up), "Downregulated:", nrow(sig_dn), "\n")

# KEGG enrichment (change organism code as needed)
kegg_up <- enrichKEGG(gene = sig_up$GeneID, organism = "hsa",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kegg_dn <- enrichKEGG(gene = sig_dn$GeneID, organism = "hsa",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Display results
cat("\n=== Upregulated KEGG pathways ===\n")
print(as.data.frame(kegg_up)[, c("Description", "Count", "p.adjust")])

cat("\n=== Downregulated KEGG pathways ===\n")
print(as.data.frame(kegg_dn)[, c("Description", "Count", "p.adjust")])