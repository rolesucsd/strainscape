#!/bin/bash

# Directory containing the input files
input_dir="../input/samples"

# File containing full paths
full_path_file="../input/fastq_files_list.txt"

# Check if the full_path.txt file exists
if [ ! -f "$full_path_file" ]; then
    echo "Error: fastq_files_list.txt file not found!"
    exit 1
fi

# Iterate over all input files in the directory
for input_file in "$input_dir"/*.txt; do
    # Get the base name of the input file (e.g., C3029.txt)
    base_name=$(basename "$input_file" .txt)

    # Output file for matched paths
    output_file="$input_dir/${base_name}_full_path.txt"

    # Ensure the input file is readable
    if [ ! -r "$input_file" ]; then
        echo "Error: Cannot read input file $input_file"
        continue
    fi

    # Process each ID in the input file
    while read -r id; do
        # Match the ID in the full_path.txt file and append to the output file
        grep "/$id" "$full_path_file" >> "$output_file"
    done < "$input_file"

    echo "Matched paths for $base_name written to $output_file"
done
