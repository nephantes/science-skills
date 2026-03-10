When analysing genomic variant data (e.g. CHIP mutations, VCF-derived tables):
- Load Excel/CSV files with readxl or tidyverse
- Filter variants by quality (PASS, DP≥10, GQ≥20 etc.) as described
- Calculate variant allele frequency (VAF) distributions per group
- Use Fisher's exact test or chi-square for categorical comparisons
- Use logistic regression (glm, family=binomial) for binary outcomes
- Use rstatix or base R for group comparisons (t-test, ANOVA, Kruskal-Wallis)
- Always report effect sizes alongside p-values