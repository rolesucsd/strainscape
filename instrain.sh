#!/bin/bash -l

#SBATCH --job-name=inStrain_profile
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

source /home/roles/anaconda3/bin/activate
conda activate sepsis

# Define the species
SPECIES="kpneumoniae"

# Define the path to the genome files
GENOME_PATH="/panfs/roles/Sepsis/${SPECIES}/drep3/dereplicated_genomes/${SPECIES}_dereplicated_genome.fna"
#GENES="/panfs/roles/Sepsis/${SPECIES}/drep/dereplicated_genomes/${SPECIES}_dereplicated_genome_cds.fna"
STB_PATH="/panfs/roles/Sepsis/${SPECIES}/drep3/dereplicated_genomes/${SPECIES}.stb"

# Output directory
OUTPUT="/panfs/roles/Sepsis/${SPECIES}/instrain"
mkdir -p "${OUTPUT}"

# Temporary directory for BAM files
TEMP_DIR="${OUTPUT}/temp_bam_files"
mkdir -p "${TEMP_DIR}"

# Calculate the file index for the current SLURM_ARRAY_TASK_ID
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Path to the file containing the full paths
FULL_PATH_FILE="/panfs/roles/Sepsis/instrain_output/full_path_list.txt"

# Read the full paths into an array
#2783 total files
mapfile -t BAM_FILES < "$FULL_PATH_FILE"
SUBSET_BAM_FILES=("${BAM_FILES[@]:2001:1000}")  # Subset to include files from index 2330 to 3330

TOTAL_FILES=${#SUBSET_BAM_FILES[@]}

if [ $FILE_INDEX -lt $TOTAL_FILES ]; then
    BAM_FILE=${SUBSET_BAM_FILES[$FILE_INDEX]}
    METAGENOMIC_FILE_NAME=$(basename "${BAM_FILE}" .bam)
    TEMP_BAM_FILE="${TEMP_DIR}/${METAGENOMIC_FILE_NAME}.bam"
    OUTPUT_DIR="${OUTPUT}/${METAGENOMIC_FILE_NAME}-vs-${SPECIES}"
    
    # Check if OUTPUT_DIR already exists
    if [ -d "$OUTPUT_DIR" ]; then
        echo "Output directory $OUTPUT_DIR already exists. Skipping..."
        exit 0  # Exit the script if OUTPUT_DIR exists
    fi
    
    # Copy BAM file to the temp directory
    cp "$BAM_FILE" "$TEMP_BAM_FILE"
    
    echo "Executing inStrain"
    # Execute the inStrain profile command
#    inStrain profile "$TEMP_BAM_FILE" "$GENOME_PATH" -o "$OUTPUT_DIR" -p 32 -g "$GENES" -s "$STB_PATH" --database_mode -d --pairing_filter all_reads
    inStrain profile "$TEMP_BAM_FILE" "$GENOME_PATH" -o "$OUTPUT_DIR" -p 32 -s "$STB_PATH" --database_mode -d --pairing_filter all_reads
else
    echo "Index ${FILE_INDEX} is out of range for available BAM files."
fi


#BAM_FILES=($(find /pscratch/dtmcdonald/sepsis-mappings-offtarget-bam/15004/ -type d -name "blongum_dereplicated_genome" | xargs -I {} find {} -type f -name "*.gz.bam"))

#python parse_stb.py --reverse -f /panfs/roles/Sepsis/blongum/drep/dereplicated_genomes/*_genomic.fna -o /panfs/roles/Sepsis/blongum/drep/dereplicated_genomes/blongum.stb

#cat /panfs/roles/Sepsis/blongum/drep/dereplicated_genomes/blongum_dereplicated_genome.fna | gzip > /projects/panpiper_pangenomes/depreplicated-genomes-for-sepsis/blongum_dereplicated_genome.fna.gz

#inStrain parse_annotations -i N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/ -o genes_output_v1 -a /panfs/roles/Sepsis/blongum/drep/dereplicated_genomes/blongum_dereplicated_genome_cds.fna

#inStrain compare -i blongum/SAM/* -o inStrain_blongum -p 6 -s /panfs/roles/Sepsis/blongum/drep/dereplicated_genomes/ --database_mode


#cat bifidobacterium/drep/data/prodigal/GCA_016415705.1_ASM1641570v1_genomic.fna.fna bifidobacterium/drep/data/prodigal/GCA_003285165.1_ASM328516v1_genomic.fna.fna  > bifidobacterium/drep/dereplicated_genomes/bifidobacterium_dereplicated_genome_cds.fna


#cat /panfs/roles/Sepsis/kpneumoniae/drep3/dereplicated_genomes/kpneumoniae_dereplicated_genome.fna | gzip > /projects/panpiper_pangenomes/depreplicated-genomes-for-sepsis/kpneumoniae_dereplicated_genome.fna.gz

#cat /panfs/roles/Sepsis/kpneumoniae/drep3/dereplicated_genomes/*.fna > /panfs/roles/Sepsis/kpneumoniae/drep3/dereplicated_genomes/kpneumoniae_dereplicated_genome.fna
