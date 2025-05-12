#!/usr/bin/env bash

# Check if input and output files are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 input_fasta output_fasta"
    exit 1
fi

input_file="$1"
output_file="$2"

awk -F'#' '
BEGIN {
    # Field separator is '#'
}
{
    if ($0 ~ /^>/) {
        # This is a header line
        # $1 is everything before the first '#'
        # Trim trailing spaces just in case
        gsub(/[[:space:]]+$/, "", $1)
        print $1
    } else {
        # Not a header line, print as is
        print $0
    }
}' "$input_file" > "$output_file"
