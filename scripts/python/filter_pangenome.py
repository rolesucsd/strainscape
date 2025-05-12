#!/bin/bash

source /home/roles/anaconda3/bin/activate

# Path to your files
PANGENOME_FA="pan_genome_reference.fa"
PANGENOME_STB="pan_genome_reference.stb"
OUTPUT_FA="pan_genome_reference_subset.fa"

# Extract gene names from the first column of pangenome.stb
awk '{print $1}' $PANGENOME_STB > gene_names.txt

# Use seqtk to subset the pangenome.fa file
seqtk subseq $PANGENOME_FA gene_names.txt > $OUTPUT_FA

# Clean up
rm gene_names.txt

echo "Subsetted pangenome.fa saved to $OUTPUT_FA"
