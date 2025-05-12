#!/bin/bash

# Directory where the genomes are stored
base_dir="/projects/panpiper_pangenomes/panpiper/BF/Assembly"

# Output files
output_file="drep.fna"
stb_file="drep.stb"

# Clear the output files if they already exist
> "$output_file"
> "$stb_file"

# Read each line from genomes.csv and process it
python3 - <<END
import os

base_dir = "${base_dir}"
output_file = "${output_file}"
stb_file = "${stb_file}"

with open('/panfs/roles/iHMP/drep/drep_selected_genomes_name.csv', 'r') as file:
    for line in file:
        sample = line.strip()
        fna_file = f"{base_dir}/{sample}/{sample}.fna"
        if os.path.isfile(fna_file):
            with open(fna_file, 'r') as fna:
                content = fna.read()
                headers = [line[1:].split()[0] for line in content.split('\n') if line.startswith('>')]
                
                # Rename headers
                renamed_content = ""
                for entry in content.split('>'):
                    if entry.strip():
                        header, *seq_lines = entry.split('\n')
                        new_header = f">{header.split()[0]}_{sample}"
                        renamed_content += new_header + "\n" + "\n".join(seq_lines) + "\n"

            with open(output_file, 'a') as out:
                out.write(renamed_content)
                
            with open(stb_file, 'a') as stb:
                for header in headers:
                    stb.write(f"{header}_{sample}\t{sample}.fna\n")

            # Save the renamed content back to the original file
            with open(fna_file, 'w') as fna:
                fna.write(renamed_content)
        else:
            print(f"File not found: {fna_file}")
END

echo "Concatenation complete. Output files: $output_file and $stb_file"
