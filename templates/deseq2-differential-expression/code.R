library(DESeq2)
library(tibble)

counts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("sample_metadata.csv")
meta$condition <- factor(meta$condition)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = meta,
  design    = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.01)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
write.csv(res_df, "de_results.csv", row.names = FALSE)
cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.01, na.rm = TRUE), "\n")