import subprocess
import os

ref = "reference.fna"
r1 = "sample_R1.fastq"
r2 = "sample_R2.fastq"

# Step 1: Trimmomatic
subprocess.run([
    "trimmomatic", "PE", r1, r2,
    "trimmed_1P.fq", "trimmed_1U.fq",
    "trimmed_2P.fq", "trimmed_2U.fq",
    "ILLUMINA:TruSeq3-PE.fa:2:30:10",
    "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
], check=True)

# Step 2: BWA alignment
subprocess.run(["bwa", "index", ref], check=True)
subprocess.run(
    f"bwa mem -R '@RG\\tID:sample\\tSM:sample\\tPL:ILLUMINA' {ref} trimmed_1P.fq trimmed_2P.fq | samtools sort -o aligned.bam",
    shell=True, check=True
)
subprocess.run(["samtools", "index", "aligned.bam"], check=True)

# Step 3: Coverage
result = subprocess.run(["samtools", "depth", "-a", "aligned.bam"], capture_output=True, text=True)
depths = [int(l.split("\t")[2]) for l in result.stdout.strip().split("\n") if l]
avg_coverage = sum(depths) / len(depths)
print(f"Average coverage: {avg_coverage:.4f}")