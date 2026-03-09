import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts = pd.read_csv("counts.csv", index_col=0)
meta = pd.read_csv("sample_metadata.csv", index_col=0)

# Ensure integer counts, genes as rows
counts = counts.round().astype(int)

# Align sample order
common = counts.columns.intersection(meta.index)
counts = counts[common]
meta = meta.loc[common]

dds = DeseqDataSet(
    counts=counts.T,  # pydeseq2 expects samples as rows
    metadata=meta,
    design="~condition",
)
dds.deseq2()

# Extract results — derive contrast from condition levels
levels = sorted(meta["condition"].unique())
contrast = ["condition", levels[0], levels[-1]]
stat_res = DeseqStats(dds, contrast=contrast)
stat_res.summary()
res_df = stat_res.results_df.sort_values("padj")

res_df.to_csv("de_results.csv")
sig = (res_df["padj"] < 0.05).sum()
print(f"Significant genes (padj < 0.05): {sig}")
print(f"Total genes tested: {len(res_df)}")