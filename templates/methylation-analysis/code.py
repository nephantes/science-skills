import pandas as pd
import numpy as np
from scipy.stats import chisquare

# Load CpG data and chromosome lengths
cpgs = pd.read_csv("cpg_sites.csv")
chrom_lengths = pd.read_csv("chromosome_lengths.csv", encoding="utf-8-sig")
print(f"CpG sites: {len(cpgs)}, unique positions: {cpgs.Pos.nunique()}")

# Filter by methylation level (hyper >90% or hypo <10%)
filtered = cpgs[(cpgs.MethylationPercentage > 90) | (cpgs.MethylationPercentage < 10)]
print(f"After filtering: {filtered.Pos.nunique()} unique sites")

# Count unique CpGs per chromosome
cpg_counts = filtered.groupby("Chromosome")["Pos"].nunique().reset_index()
cpg_counts.columns = ["Chromosome", "n_cpgs"]

# Merge with chromosome lengths
chrom_lengths["Chromosome"] = chrom_lengths["Chromosome"].astype(str)
cpg_counts["Chromosome"] = cpg_counts["Chromosome"].astype(str)
merged = cpg_counts.merge(chrom_lengths, on="Chromosome", how="inner")

# Compute density
merged["density"] = merged["n_cpgs"] / merged["Length"]
print(f"\nPer-chromosome density:")
print(merged.sort_values("Chromosome").to_string(index=False))

# Genome-wide average density
avg_density = merged["density"].mean()
print(f"\nGenome-wide average density: {avg_density:.4e}")

# Chi-square goodness-of-fit: is CpG distribution proportional to chromosome length?
observed = merged["n_cpgs"].values
expected_prop = merged["Length"].values / merged["Length"].values.sum()
expected = expected_prop * observed.sum()
chi2, p = chisquare(observed, f_exp=expected)
print(f"\nChi-square test: chi2={chi2:.2f}, p={p:.4e}")