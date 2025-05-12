#!/bin/bash -l

#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=50G

source /home/roles/anaconda3/bin/activate
conda activate mapping

bwa index /projects/panpiper_pangenomes/panpiper/BF2/Bacteroides/GCF_018292165.1_ASM1829216v1_genomic.fna

#bowtie2-build /projects/panpiper_pangenomes/panpiper/BF2/Bacteroides/GCF_018292165.1_ASM1829216v1_genomic.fna /projects/panpiper_pangenomes/panpiper/BF2/Bacteroides