import pandas as pd
import numpy as np
from scipy import stats

df = pd.read_csv("data.csv")
print(f"Data shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
print(f"Groups: {df['group'].unique()}")

# Example: Two-sample t-test (Welch's)
group_a = df[df['group'] == 'A']['value']
group_b = df[df['group'] == 'B']['value']
t_stat, p_val = stats.ttest_ind(group_a, group_b, equal_var=False)
print(f"Welch's t-test: t={t_stat:.4f}, p={p_val:.4e}")

# Cohen's d (effect size)
pooled_std = np.sqrt((group_a.std()**2 + group_b.std()**2) / 2)
cohens_d = (group_a.mean() - group_b.mean()) / pooled_std
print(f"Cohen's d: {cohens_d:.4f}")

# Shapiro-Wilk normality test
for grp in df['group'].unique():
    vals = df[df['group'] == grp]['value']
    stat, p = stats.shapiro(vals)
    print(f"Shapiro-Wilk ({grp}): W={stat:.4f}, p={p:.4e}")

# One-way ANOVA (if >2 groups)
groups = [g['value'].values for _, g in df.groupby('group')]
if len(groups) > 2:
    f_stat, p_anova = stats.f_oneway(*groups)
    print(f"ANOVA: F={f_stat:.4f}, p={p_anova:.4e}")