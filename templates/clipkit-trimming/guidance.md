When trimming alignments with ClipKIT:
- ClipKIT is a command-line tool: 'clipkit input.fasta' (default mode is smart-gap)
- Output file is automatically named input.fasta.clipkit
- Compare alignment length before and after: parse both with Bio.AlignIO
- Percentage reduction = (original_length - trimmed_length) / original_length × 100
- Report per-file statistics and identify which file has the highest/lowest reduction
- Ortholog ID is typically the filename prefix before "at" (e.g., "1014314" from "1014314at2759.faa.mafft")