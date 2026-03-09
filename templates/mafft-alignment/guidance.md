When running MAFFT alignment on FASTA files:
- Use subprocess to call 'mafft --auto input.faa' and capture stdout as the aligned FASTA
- Parse aligned output with Bio.AlignIO (format="fasta")
- Alignment length = number of columns in the alignment (all sequences same length after alignment)
- Gap character is "-"; count gaps across ALL sequences in ALL alignments
- Total positions = sum of (alignment_length × num_sequences) for each alignment
- Gap percentage = (total_gaps / total_positions) × 100
- Save each aligned output as inputname.mafft for downstream use
- Report per-file and aggregate statistics