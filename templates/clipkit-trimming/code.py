import subprocess
import glob
import os
from Bio import AlignIO

aligned_files = sorted(glob.glob("*.mafft") + glob.glob("*.aln") + glob.glob("*.fasta"))
print(f"Found {len(aligned_files)} alignment files to trim")

results = []
for aln_file in aligned_files:
    # Measure original length
    aln_orig = AlignIO.read(aln_file, "fasta")
    orig_len = aln_orig.get_alignment_length()

    # Run ClipKIT
    subprocess.run(["clipkit", aln_file], check=True)
    trimmed_file = aln_file + ".clipkit"

    # Measure trimmed length
    aln_trimmed = AlignIO.read(trimmed_file, "fasta")
    trim_len = aln_trimmed.get_alignment_length()

    reduction = (orig_len - trim_len) / orig_len * 100 if orig_len > 0 else 0
    # Extract ortholog ID from filename
    ortholog_id = os.path.basename(aln_file).split("at")[0]

    results.append({
        "ortholog_id": ortholog_id,
        "file": aln_file,
        "original_length": orig_len,
        "trimmed_length": trim_len,
        "reduction_pct": round(reduction, 2),
    })
    print(f"  {ortholog_id}: {orig_len} -> {trim_len} ({reduction:.1f}% reduction)")

import pandas as pd
df = pd.DataFrame(results).sort_values("reduction_pct", ascending=False)
df.to_csv("trimming_stats.csv", index=False)

top = df.iloc[0]
print(f"\nHighest reduction: ortholog {top['ortholog_id']} ({top['reduction_pct']:.1f}%)")