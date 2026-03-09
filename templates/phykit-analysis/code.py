import subprocess
import glob
import os
import numpy as np

treefiles = sorted(glob.glob("*.treefile"))
print(f"Found {len(treefiles)} tree files")

results = []
for tf in treefiles:
    ortholog_id = os.path.basename(tf).split("at")[0]
    result = subprocess.run(["phykit", "treeness", tf], capture_output=True, text=True)
    value = float(result.stdout.strip())
    results.append({"ortholog_id": ortholog_id, "file": tf, "treeness": value})
    print(f"  {ortholog_id}: treeness = {value:.4f}")

values = [r["treeness"] for r in results]
avg = np.mean(values)
median = np.median(values)
print(f"\nAverage treeness: {avg:.4f}")
print(f"Median treeness: {median:.4f}")
print(f"Average × 1000 (rounded): {round(avg * 1000)}")

import pandas as pd
pd.DataFrame(results).to_csv("phykit_treeness.csv", index=False)