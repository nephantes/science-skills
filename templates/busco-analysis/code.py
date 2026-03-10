import subprocess
import os
import glob

# Run BUSCO on a proteome
proteome = "species.faa"
lineage = "eukaryota_odb10"
out_name = "busco_output"

subprocess.run([
    "busco", "-i", proteome, "-m", "protein",
    "-l", lineage, "-o", out_name, "--force"
], check=True)

# Parse summary
summary_file = glob.glob(f"{out_name}/short_summary*.txt")[0]
with open(summary_file) as f:
    for line in f:
        if "C:" in line and "S:" in line:
            print(line.strip())

# Find single-copy orthologs
sco_dir = f"{out_name}/run_{lineage}/busco_sequences/single_copy_busco_sequences/"
sco_files = glob.glob(os.path.join(sco_dir, "*.faa"))
print(f"Single-copy orthologs found: {len(sco_files)}")

# Count total amino acids across all orthologs
total_aa = 0
for faa in sco_files:
    with open(faa) as f:
        for line in f:
            if not line.startswith(">"):
                total_aa += len(line.strip())
print(f"Total amino acids: {total_aa}")