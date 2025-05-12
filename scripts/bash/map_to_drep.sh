#!/bin/bash -l

#SBATCH --partition=short
#SBATCH --job-name=process-fastq
#SBATCH --output=process_fastq_%A_%a.out
#SBATCH --error=process_fastq_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-21%21

# Activate environment
source /home/roles/anaconda3/bin/activate
conda activate mapping

# Read the FASTQ file path from the list
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /panfs/roles/iHMP/input/H4020_full_path.txt)
SAMPLE=$(basename $FASTQ .fastq.gz)

# Define output directories
ALIGNMENT_DIR="/panfs/roles/iHMP/H4020/alignment_pangenome"

# Make sure output directories exist
mkdir -p $ALIGNMENT_DIR

# Processing steps
OUTPUT_FILE="$ALIGNMENT_DIR/${SAMPLE}.sam"

# Alignment
if [ ! -f "$OUTPUT_FILE" ]; then
    gzip -dc $FASTQ | bowtie2 -p 8 -x /panfs/roles/iHMP/drep/pan_genome_reference -U - -S $OUTPUT_FILE
fi

echo "Processing complete for $FASTQ"
