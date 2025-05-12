#!/bin/bash -l

#SBATCH --partition=short
#SBATCH --job-name=process-fastq
#SBATCH --output=process_fastq_%A_%a.out
#SBATCH --error=process_fastq_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-614%50

# Activate environment
source /home/roles/anaconda3/bin/activate
conda activate mapping

# Read the FASTQ file path from the list
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /ddn_scratch/roles/strain_analysis/iHMP/input/fastq_files_list2.txt)
SAMPLE=$(basename $FASTQ .fastq.gz)

# Define output directories
ALIGNMENT_DIR="/ddn_scratch/roles/strain_analysis/iHMP/coverages/BU/alignment_results"
COVERAGE_DIR="/ddn_scratch/roles/strain_analysis/iHMP/coverages/BU/coverage_results"

# Make sure output directories exist
mkdir -p $ALIGNMENT_DIR $COVERAGE_DIR

# Paths
SORTED_BAM="$ALIGNMENT_DIR/${SAMPLE}_sorted.bam"
COVERAGE_FILE="$COVERAGE_DIR/${SAMPLE}_coverage.txt"
DEPTH_FILE="$COVERAGE_DIR/${SAMPLE}_depth.txt"
REFERENCE="/projects/panpiper_pangenomes/panpiper/BF2/Bacteroides/GCF_018292165.1_ASM1829216v1_genomic.fna"

# Alignment and conversion directly to BAM
if [ ! -f "$SORTED_BAM" ]; then
    gzip -dc $FASTQ | bwa mem -t 8 $REFERENCE - |
    samtools view -bS - | samtools sort -o $SORTED_BAM
    samtools index $SORTED_BAM
fi

# Coverage calculation
if [ ! -f "$COVERAGE_FILE" ]; then
    samtools depth $SORTED_BAM > $COVERAGE_FILE
fi