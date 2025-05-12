#!/bin/bash -l

#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=30G

# Array of directories containing FASTQ files
DIRS=(90155 90164 90165 90462 90463 100621 100620)

# Activate environment
source /home/roles/anaconda3/bin/activate
conda activate bwa

# Directory for the specific array job
FASTQ_DIR="/qmounts/qiita_data/per_sample_FASTQ/${DIRS[0]}"
ALIGNMENT_DIR="alignment_results"
COVERAGE_DIR="coverage_results"

# Process each FASTQ file in the directory
for FASTQ in $FASTQ_DIR/*.fastq.gz
do
    SAMPLE=$(basename $FASTQ .fastq.gz)
    SAM_FILE="$ALIGNMENT_DIR/${SAMPLE}.sam"
    
    # Check if the alignment file already exists
    if [ ! -f "$SAM_FILE" ]; then
        gzip -dc $FASTQ | bwa mem -t 8 /projects/panpiper_pangenomes/panpiper/BF/reference/9343.fna - > $SAM_FILE
    fi
done

# Deactivate BWA environment
conda deactivate
conda activate concoct

# Further processing: coverage results
for SAM_FILE in $ALIGNMENT_DIR/*.sam
do
    SAMPLE=$(basename $SAM_FILE .sam)
    BAM_FILE="$ALIGNMENT_DIR/${SAMPLE}.bam"
    SORTED_BAM="$ALIGNMENT_DIR/${SAMPLE}_sorted.bam"
    COVERAGE_FILE="$COVERAGE_DIR/${SAMPLE}_coverage.txt"
    
    # Convert SAM to BAM, sort, index, and calculate coverage if not already done
    if [ ! -f "$SORTED_BAM" ]; then
        samtools view -bS $SAM_FILE > $BAM_FILE
        samtools sort $BAM_FILE -o $SORTED_BAM
        samtools index $SORTED_BAM
    fi
    
    # Check if the coverage file already exists
    if [ ! -f "$COVERAGE_FILE" ]; then
        samtools depth -a $SORTED_BAM > $COVERAGE_FILE
    fi
done

echo "Processing complete for directory ${DIRS[$SLURM_ARRAY_TASK_ID]}"
