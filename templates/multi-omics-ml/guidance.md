When performing multi-omics integration or ML analysis on biological data:
- Load gene expression, methylation, and clinical data separately
- Handle missing data: count complete cases across modalities
- Compute distribution statistics (skewness, kurtosis) to characterise data shape
- Use scipy.stats.skew or pandas .skew() for skewness assessment
- For classification: use sklearn (LogisticRegression, RandomForest) or xgboost
- For survival prediction: use lifelines or sklearn with Cox regression features
- Always split data before fitting (train_test_split with random_state for reproducibility)
- Report confusion matrix, AUC, or other appropriate metrics