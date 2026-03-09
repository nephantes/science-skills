import subprocess
import glob
import os
from Bio import AlignIO
from io import StringIO

fasta_files = sorted(glob.glob("*.faa"))
print(f"Found {len(fasta_files)} FASTA files to align")

total_gaps = 0
total_positions = 0
results = []

for fasta in fasta_files:
    # Run MAFFT
    result = subprocess.run(
        ["mafft", "--auto", fasta],
        capture_output=True, text=True
    )
    aligned_text = result.stdout

    # Save aligned output
    out_name = fasta + ".mafft"
    with open(out_name, "w") as f:
        f.write(aligned_text)

    # Parse alignment
    aln = AlignIO.read(StringIO(aligned_text), "fasta")
    aln_len = aln.get_alignment_length()
    n_seqs = len(aln)
    gaps = sum(str(rec.seq).count("-") for rec in aln)
    positions = aln_len * n_seqs

    total_gaps += gaps
    total_positions += positions

    pct = (gaps / positions * 100) if positions > 0 else 0
    results.append({"file": fasta, "seqs": n_seqs, "aln_len": aln_len, "gaps": gaps, "gap_pct": round(pct, 2)})
    print(f"  {fasta}: {n_seqs} seqs, {aln_len} cols, {gaps} gaps ({pct:.2f}%)")

overall_pct = (total_gaps / total_positions * 100) if total_positions > 0 else 0
print(f"\nOverall: {total_gaps} gaps / {total_positions} positions = {overall_pct:.1f}%")

import pandas as pd
pd.DataFrame(results).to_csv("alignment_stats.csv", index=False)