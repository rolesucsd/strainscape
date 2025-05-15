#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=25GB
#SBATCH --cpus-per-task=8

source /home/roles/anaconda3/bin/activate

# List of sample directories
SAMPLES=$(find /ddn_scratch/roles/strain_analysis/iHMP/instrain_output -maxdepth 1 -type d)

# Iterate over each sample directory
while IFS= read -r SAMPLE; do
    echo "Processing sample: $SAMPLE"
    python -u iHMP_scripts/iHMP_main.py --file_dir "$SAMPLE"
done <<< "$SAMPLES"