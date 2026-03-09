When performing PCA, always: (1) transpose so samples are rows,
(2) remove zero-variance features, (3) cap n_components to min(n_samples-1, n_features),
(4) color by sample groups from sample_metadata.csv.