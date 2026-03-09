When performing PCA in R, always: (1) transpose so samples are rows,
(2) remove zero-variance features with check via apply(X, 2, var),
(3) use prcomp with center=TRUE and scale.=TRUE,
(4) color by sample groups from sample_metadata.csv.