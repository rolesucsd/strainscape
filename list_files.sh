#!/bin/bash

BASEPATH="/pscratch/dtmcdonald/sepsis-mappings-offtarget-bam/15004"

IDENTIFIERS_FILE="metadata_sample_id.txt"

OUTPUT_FILE="instrain_output/full_path_list.txt"

> "$OUTPUT_FILE"

while IFS= read -r identifier; do
    find "$BASEPATH" -type f -path "*/blongum_dereplicated_genome/${identifier}_*_L001_R1_001.trimmed.fastq.gz.bam" >> "$OUTPUT_FILE"
done < "$IDENTIFIERS_FILE"

echo "Found file paths written to $OUTPUT_FILE"
