suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
library(tibble)
library(ggplot2)
library(reshape2)
library(tidyr)

# Load counts — first column is gene IDs
counts_raw <- read.csv("counts.csv", check.names = FALSE)
ids <- make.unique(as.character(counts_raw[[1]]), sep = "_")
counts <- counts_raw[, -1, drop = FALSE]
rownames(counts) <- ids

# Load sample metadata — must have a 'sample' column and at least one group column
meta <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)
group_col <- setdiff(colnames(meta), "sample")[1]
meta$group <- factor(meta[[group_col]])
rownames(meta) <- meta$sample

# Align samples
common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
cat("Counts shape:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Groups:", paste(levels(meta$group), collapse = ", "), "\n")

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")
dds <- DESeq(dds)

# All pairwise comparisons
groups <- levels(meta$group)
comparisons <- combn(groups, 2, simplify = FALSE)
cat("Number of pairwise comparisons:", length(comparisons), "\n")
dir.create("de_results", showWarnings = FALSE)

summary_df <- data.frame(
  comparison = character(), group1 = character(), group2 = character(),
  total_genes = integer(), sig_genes_padj05 = integer(),
  up_in_group2 = integer(), down_in_group2 = integer(),
  stringsAsFactors = FALSE
)

for (comp in comparisons) {
  g1 <- comp[1]; g2 <- comp[2]
  comp_name <- paste(g2, "vs", g1, sep = "_")
  res <- results(dds, contrast = c("group", g2, g1))
  res <- res[order(res$pvalue), ]
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  write.table(res_df, file = paste0("de_results/", comp_name, ".tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  up  <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  dn  <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  summary_df <- rbind(summary_df, data.frame(
    comparison = comp_name, group1 = g1, group2 = g2,
    total_genes = nrow(res_df), sig_genes_padj05 = sig,
    up_in_group2 = up, down_in_group2 = dn
  ))
  cat(sprintf("Completed: %s - %d significant genes\n", comp_name, sig))
}

write.table(summary_df, "de_results/summary_all_comparisons.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)
print(summary_df)

# Visualization — always create before any downstream plots
if (nrow(summary_df) > 1) {
  # 1. Pairwise DEG count heatmap
  n <- length(groups)
  deg_matrix <- matrix(0L, nrow = n, ncol = n)
  rownames(deg_matrix) <- groups
  colnames(deg_matrix) <- groups
  for (i in seq_len(nrow(summary_df))) {
    g1 <- summary_df$group1[i]; g2 <- summary_df$group2[i]
    deg_matrix[g1, g2] <- summary_df$sig_genes_padj05[i]
    deg_matrix[g2, g1] <- summary_df$sig_genes_padj05[i]
  }
  deg_long <- melt(deg_matrix)
  colnames(deg_long) <- c("Group1", "Group2", "DEGs")
  p_heat <- ggplot(deg_long, aes(x = Group1, y = Group2, fill = DEGs)) +
    geom_tile(color = "white") +
    geom_text(aes(label = DEGs), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "#1f77b4",
                        name = "Significant\nDEGs\n(padj<0.05)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    labs(title = "Pairwise Differential Expression Analysis\n(Number of Significant DEGs, padj < 0.05)") +
    coord_fixed()
  ggsave("de_results/deg_heatmap.pdf", p_heat, width = 10, height = 8)
  ggsave("de_results/deg_heatmap.png", p_heat, width = 10, height = 8, dpi = 300)

  # 2. Up/Down direction barplot
  comp_order <- summary_df %>% arrange(desc(sig_genes_padj05)) %>% pull(comparison)
  summary_long <- summary_df %>%
    select(comparison, up_in_group2, down_in_group2) %>%
    pivot_longer(cols = c(up_in_group2, down_in_group2),
                 names_to = "direction", values_to = "count") %>%
    mutate(direction = ifelse(direction == "up_in_group2", "Up-regulated", "Down-regulated"),
           count = ifelse(direction == "Down-regulated", -count, count),
           comparison = factor(comparison, levels = rev(comp_order)))
  bar_height <- max(6, nrow(summary_df) * 0.4 + 2)
  p_bar <- ggplot(summary_long, aes(x = comparison, y = count, fill = direction)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("Up-regulated" = "#e74c3c", "Down-regulated" = "#3498db")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8), legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    labs(title = "Differentially Expressed Genes by Comparison",
         x = "Comparison", y = "Number of DEGs (padj < 0.05)", fill = "Direction") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
  ggsave("de_results/deg_barplot.pdf", p_bar, width = 10, height = bar_height)
  ggsave("de_results/deg_barplot.png", p_bar, width = 10, height = bar_height, dpi = 300)
} else {
  # Single comparison — volcano plot
  res_df <- read.table(paste0("de_results/", summary_df$comparison[1], ".tsv"),
    header = TRUE, sep = "\t")
  res_df$significant <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  p_volc <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#e74c3c")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = paste("Volcano:", summary_df$comparison[1]),
         x = "log2 Fold Change", y = "-log10(p-value)", color = "Significant") +
    theme_minimal(base_size = 12)
  ggsave("de_results/volcano_plot.png", p_volc, width = 8, height = 6, dpi = 300)
}
