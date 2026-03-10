When performing BUSCO analysis:
- Run BUSCO via subprocess: busco -i <proteome.faa> -m protein -l <lineage> -o <output_name>
- For genome mode: busco -i <genome.fna> -m genome -l <lineage>
- Parse the short_summary*.txt file for completeness stats (Complete, Single, Duplicated, Fragmented, Missing)
- To find shared single-copy orthologs across species, compare run_*/busco_sequences/single_copy_busco_sequences/
- Count amino acids in ortholog FASTA files using BioPython or simple parsing
- Common lineages: eukaryota_odb10, metazoa_odb10, fungi_odb10, vertebrata_odb10
- If lineage dataset is provided as .tar.gz, extract it first and use --offline flag