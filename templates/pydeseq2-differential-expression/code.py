import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os

counts = pd.read_csv("counts.csv", index_col=0)
meta = pd.read_csv("sample_metadata.csv", index_col=0)

# Ensure integer counts; genes as rows, samples as columns
counts = counts.round().astype(int)

# Use first metadata column as group variable
group_col = meta.columns[0]
meta[group_col] = meta[group_col].astype(str)
print(f"Group column: '{group_col}', groups: {sorted(meta[group_col].unique())}")

# Align samples
common = counts.columns.intersection(meta.index)
counts = counts[common]
meta = meta.loc[common]
print(f"Counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

# Run DESeq2 via pydeseq2
dds = DeseqDataSet(
    counts=counts.T,  # pydeseq2 expects samples as rows
    metadata=meta,
    design_factors=group_col,
)
dds.deseq2()

groups = sorted(meta[group_col].unique())
comps = list(combinations(groups, 2))
print(f"Number of pairwise comparisons: {len(comps)}")
os.makedirs("de_results", exist_ok=True)

summary_rows = []
for g1, g2 in comps:
    comp_name = f"{g2}_vs_{g1}"
    stat_res = DeseqStats(dds, contrast=[group_col, g2, g1])
    stat_res.summary()
    res_df = stat_res.results_df.copy()
    res_df.index.name = "gene"
    res_df = res_df.reset_index().sort_values("padj")
    res_df.to_csv(f"de_results/{comp_name}.csv", index=False)
    sig = int((res_df["padj"] < 0.05).sum())
    up  = int(((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] > 0)).sum())
    dn  = int(((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] < 0)).sum())
    summary_rows.append(dict(comparison=comp_name, group1=g1, group2=g2,
        total_genes=len(res_df), sig_genes_padj05=sig,
        up_in_group2=up, down_in_group2=dn))
    print(f"Completed: {comp_name} - {sig} significant genes")

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("de_results/summary_all_comparisons.csv", index=False)
print(summary_df.to_string(index=False))

# Visualization — always create before any downstream plots
if len(summary_df) > 1:
    # 1. Pairwise DEG count heatmap
    deg_matrix = pd.DataFrame(0, index=groups, columns=groups, dtype=int)
    for _, row in summary_df.iterrows():
        deg_matrix.loc[row["group1"], row["group2"]] = row["sig_genes_padj05"]
        deg_matrix.loc[row["group2"], row["group1"]] = row["sig_genes_padj05"]

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(deg_matrix, annot=True, fmt="d", cmap="Blues",
                linewidths=0.5, ax=ax,
                cbar_kws={"label": "Significant DEGs (padj < 0.05)"})
    ax.set_title("Pairwise Differential Expression Analysis\n(Number of Significant DEGs, padj < 0.05)",
                 fontsize=14, fontweight="bold", pad=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()
    plt.savefig("de_results/deg_heatmap.pdf", bbox_inches="tight")
    plt.savefig("de_results/deg_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. Up/Down direction barplot
    comp_order = summary_df.sort_values("sig_genes_padj05", ascending=False)["comparison"].tolist()
    bar_data = []
    for _, row in summary_df.iterrows():
        bar_data.append({"comparison": row["comparison"], "count": row["up_in_group2"], "direction": "Up-regulated"})
        bar_data.append({"comparison": row["comparison"], "count": -row["down_in_group2"], "direction": "Down-regulated"})
    bar_df = pd.DataFrame(bar_data)
    bar_df["comparison"] = pd.Categorical(bar_df["comparison"],
                                           categories=list(reversed(comp_order)), ordered=True)
    bar_height = max(6, len(summary_df) * 0.4 + 2)
    fig, ax = plt.subplots(figsize=(10, bar_height))
    colors = {"Up-regulated": "#e74c3c", "Down-regulated": "#3498db"}
    for direction, grp in bar_df.groupby("direction", observed=True):
        ax.barh(grp["comparison"].astype(str), grp["count"], color=colors[direction], label=direction)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Number of DEGs (padj < 0.05)", fontsize=10)
    ax.set_ylabel("Comparison", fontsize=10)
    ax.set_title("Differentially Expressed Genes by Comparison", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right")
    ax.tick_params(axis="y", labelsize=8)
    plt.tight_layout()
    plt.savefig("de_results/deg_barplot.pdf", bbox_inches="tight")
    plt.savefig("de_results/deg_barplot.png", dpi=300, bbox_inches="tight")
    plt.close()
else:
    # Single comparison — volcano plot
    comp_file = f"de_results/{summary_df['comparison'].iloc[0]}.csv"
    res_df = pd.read_csv(comp_file)
    sig = (res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() > 1)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(res_df.loc[~sig, "log2FoldChange"], -np.log10(res_df.loc[~sig, "pvalue"]),
               c="grey", alpha=0.5, s=10, label="Not significant")
    ax.scatter(res_df.loc[sig, "log2FoldChange"], -np.log10(res_df.loc[sig, "pvalue"]),
               c="#e74c3c", alpha=0.7, s=12, label="Significant")
    ax.axvline(-1, linestyle="--", color="grey", linewidth=0.8)
    ax.axvline(1, linestyle="--", color="grey", linewidth=0.8)
    ax.axhline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Volcano: {summary_df['comparison'].iloc[0]}")
    ax.legend()
    plt.tight_layout()
    plt.savefig("de_results/volcano_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
