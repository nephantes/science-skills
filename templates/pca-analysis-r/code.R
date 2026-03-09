library(ggplot2)

expr <- read.csv("expression_matrix.csv", row.names = 1)
meta <- read.csv("sample_metadata.csv")

# Transpose: samples as rows, genes as columns
X <- t(expr)

# Remove zero-variance genes
gene_var <- apply(X, 2, var)
X <- X[, gene_var > 0]
stopifnot(all(is.finite(X)))

# PCA
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
var_explained <- summary(pca_res)$importance[2, 1:2] * 100

# Build plot data
pca_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  condition = meta$condition
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal()
ggsave("pca_plot.png", p, width = 8, height = 6, dpi = 100)