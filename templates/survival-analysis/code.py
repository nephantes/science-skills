import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

df = pd.read_csv("clinical_data.csv")
# Expect columns: time, event (1=event, 0=censored), group

kmf = KaplanMeierFitter()
fig, ax = plt.subplots(figsize=(10, 6))

for grp in df['group'].unique():
    mask = df['group'] == grp
    kmf.fit(df.loc[mask, 'time'], df.loc[mask, 'event'], label=grp)
    kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_xlabel("Time")
ax.set_ylabel("Survival Probability")
ax.set_title("Kaplan-Meier Survival Curves")
plt.savefig("km_curves.png", dpi=100, bbox_inches="tight")
plt.close()

# Log-rank test between groups
groups = df['group'].unique()
if len(groups) == 2:
    g1 = df[df['group'] == groups[0]]
    g2 = df[df['group'] == groups[1]]
    lr = logrank_test(g1['time'], g2['time'], g1['event'], g2['event'])
    print(f"Log-rank test p-value: {lr.p_value:.4e}")