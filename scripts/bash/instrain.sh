#!/bin/bash -l

#SBATCH --job-name=inStrain_profile
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50G
#SBATCH --array=1-21%21  # Modify this based on the total number of BAM files

source /home/roles/anaconda3/bin/activate
conda activate sepsis

# Define the path to the genome files
GENOME_PATH="/panfs/roles/iHMP/drep/pan_genome_reference_subset.fa"
STB_PATH="/panfs/roles/iHMP/drep/pan_genome_reference.stb"

# Output directory
OUTPUT="/panfs/roles/iHMP/H4020/instrain_pangenome"
mkdir -p "${OUTPUT}"

# Get the list of BAM files
BAM_FILES=(/panfs/roles/iHMP/H4020/alignment_pangenome/*.sam)

# Calculate the file index for the current SLURM_ARRAY_TASK_ID
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

TOTAL_FILES=${#BAM_FILES[@]}

if [ $FILE_INDEX -lt $TOTAL_FILES ]; then
    BAM_FILE=${BAM_FILES[$FILE_INDEX]}
    METAGENOMIC_FILE_NAME=$(basename "${BAM_FILE}" .sam)
    OUTPUT_DIR="${OUTPUT}/${METAGENOMIC_FILE_NAME}"
    
    # Check if OUTPUT_DIR already exists
    if [ -d "$OUTPUT_DIR" ]; then
        echo "Output directory $OUTPUT_DIR already exists. Skipping..."
        exit 0  # Exit the script if OUTPUT_DIR exists
    fi
    
    echo "Executing inStrain"
    # Execute the inStrain profile command
    inStrain profile "$BAM_FILE" "$GENOME_PATH" -s "$STB_PATH" -o "$OUTPUT_DIR" -c 2 -p 32 --database_mode -d --pairing_filter all_reads
else
    echo "Index ${FILE_INDEX} is out of range for available BAM files."
fi
