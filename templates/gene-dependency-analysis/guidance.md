When analysing gene dependency/essentiality data (e.g. DepMap):
- Gene effect scores: negative = essential, 0 = non-essential
- Expression data: typically log2(TPM+1), often right-skewed
- For each gene, compute Spearman correlation between expression and dependency
- Apply multiple testing correction (BH/FDR) on p-values
- Identify genes where higher expression correlates with greater dependency (negative correlation)
- Report distribution characteristics (skewness) of expression data
- Use scipy.stats.spearmanr, statsmodels.stats.multitest.multipletests