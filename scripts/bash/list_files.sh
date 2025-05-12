#!/bin/bash

# Base directory containing all subdirectories with FASTQ files
ROOT_DIR="/qmounts/qiita_data/per_sample_FASTQ"

# Array of directories to include
DIRS=(90155 90164 90165 90462 90463 100621 100620)

# File to store the list of FASTQ files
FASTQ_LIST="fastq_files_list.txt"

# Clear or create the FASTQ list file
> $FASTQ_LIST

# Loop through the specified directories and add all FASTQ.gz files to the list
for DIR in "${DIRS[@]}"; do
    find "$ROOT_DIR/$DIR" -name "*.fastq.gz" >> $FASTQ_LIST
done
