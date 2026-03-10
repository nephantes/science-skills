suppressPackageStartupMessages(library(DESeq2))
library(tibble)

counts_raw <- read.csv("counts.csv", check.names = FALSE)
ids <- make.unique(as.character(counts_raw[[1]]), sep = "_")
counts <- counts_raw[, -1, drop = FALSE]
rownames(counts) <- ids

meta <- read.csv("sample_metadata.csv")
meta$condition <- factor(meta$condition)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = meta,
  design    = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
write.csv(res_df, "de_results.csv", row.names = FALSE)
cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")