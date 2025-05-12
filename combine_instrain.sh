#!/bin/bash

# Define directories and output files
INSTRAIN_OUTPUT_DIR="/panfs/roles/iHMP/M2069/instrain_9343"  # Replace with the actual path
COMBINED_OUTPUT_SCAFFOLD="/panfs/roles/iHMP/M2069/combined_scaffold_info_9343.tsv"  # Replace with the actual path
COMBINED_OUTPUT_SNV="/panfs/roles/iHMP/M2069/combined_SNV_info_9343.tsv"  # Replace with the actual path

# Get all directory names in the INSTRAIN_OUTPUT_DIR
samples=($(ls -d $INSTRAIN_OUTPUT_DIR/*/ | xargs -n 1 basename))

# Convert the samples array to a comma-separated string
samples_str=$(printf ",\"%s\"" "${samples[@]}")
samples_str="[${samples_str:1}]"

# Combine scaffold info
echo "Combining scaffold info..."
python - <<EOF
import pandas as pd
import os

# Define directories and output files
INSTRAIN_OUTPUT_DIR = "$INSTRAIN_OUTPUT_DIR"
COMBINED_OUTPUT_SCAFFOLD = "$COMBINED_OUTPUT_SCAFFOLD"

# Define samples
samples = $samples_str

# Prepare to collect data
dataframes = []
for sample in samples:
    file_path = os.path.join(INSTRAIN_OUTPUT_DIR, f"{sample}/output/{sample}_scaffold_info.tsv")
    if os.path.exists(file_path):
        try:
            df = pd.read_csv(file_path, sep='\t')
            if not df.empty:
                df['External.ID'] = sample
                dataframes.append(df)
            else:
                print(f"Skipping empty file: {file_path}")
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {file_path}")

# Combine data if available
if dataframes:
    combined_df = pd.concat(dataframes, ignore_index=True)
    combined_df.to_csv(COMBINED_OUTPUT_SCAFFOLD, sep='\t', index=False)
    print(f"Scaffold info combined and written to {COMBINED_OUTPUT_SCAFFOLD}.")
else:
    print(f"No valid scaffold files found. No output created.")
EOF

# Combine SNV info
echo "Combining SNV info..."
python - <<EOF
import pandas as pd
import os

# Define directories and output files
INSTRAIN_OUTPUT_DIR = "$INSTRAIN_OUTPUT_DIR"
COMBINED_OUTPUT_SNV = "$COMBINED_OUTPUT_SNV"

# Define samples
samples = $samples_str

# Prepare to collect data
dataframes = []
for sample in samples:
    file_path = os.path.join(INSTRAIN_OUTPUT_DIR, f"{sample}/output/{sample}_SNVs.tsv")
    if os.path.exists(file_path):
        try:
            df = pd.read_csv(file_path, sep='\t')
            if not df.empty:
                df['External.ID'] = sample
                dataframes.append(df)
            else:
                print(f"Skipping empty file: {file_path}")
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {file_path}")

# Combine data if available
if dataframes:
    combined_df = pd.concat(dataframes, ignore_index=True)
    combined_df.to_csv(COMBINED_OUTPUT_SNV, sep='\t', index=False)
    print(f"SNV info combined and written to {COMBINED_OUTPUT_SNV}.")
else:
    print(f"No valid SNV files found. No output created.")
EOF
