#!/bin/bash

# Define the main directory and output file
MAIN_DIR="/ddn_scratch/roles/strain_analysis/iHMP/reference/bakta"
OUTPUT_FILE="${MAIN_DIR}/wolr2_reference.tsv"

# Define the desired header
HEADER="Sequence Id	Type	Start	Stop	Strand	Locus Tag	Gene	Product	DbXrefs"

# Check if the main directory exists
if [ ! -d "$MAIN_DIR" ]; then
    echo "Error: Directory $MAIN_DIR does not exist."
    exit 1
fi

# Write the header to the output file
echo -e "$HEADER" > "$OUTPUT_FILE"

# Find all {sample}.tsv files in subdirectories and append their data (excluding headers)
find "$MAIN_DIR" -type f -name "*.tsv" | while read -r file; do
    # Skip files that end with .hypotheticals.tsv or .inference.tsv
    if [[ "$file" =~ \.hypotheticals\.tsv$ || "$file" =~ \.inference\.tsv$ ]]; then
        continue
    fi
    
    # Check if the file has the expected Bakta header and skip those lines
    if grep -q "# Annotated with Bakta" "$file"; then
        # Append content after skipping the header lines (first 6 lines)
        tail -n +7 "$file" >> "$OUTPUT_FILE"
    else
        echo "Warning: $file does not appear to have the expected Bakta header. Skipping."
    fi
done

# Check if the output file was created successfully
if [ -f "$OUTPUT_FILE" ]; then
    echo "Successfully combined .tsv files into $OUTPUT_FILE"
else
    echo "Error: Failed to create $OUTPUT_FILE"
    exit 1
fi