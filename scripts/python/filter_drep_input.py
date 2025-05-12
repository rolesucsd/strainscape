# Define the path to the files
file_to_filter = '/panfs/roles/iHMP/input/full_path_list.txt'
values_file = '/panfs/roles/iHMP/input/M2069.txt'
output_file = '/panfs/roles/iHMP/input/M2069_bam.txt'

# Read the values to filter by into a set (for faster lookup)
with open(values_file, 'r') as vf:
    values_to_filter = set(line.strip() for line in vf if line.strip())

# Open the output file to write the filtered lines
with open(output_file, 'w') as of:
    # Open the file to be filtered
    with open(file_to_filter, 'r') as df:
        # Iterate over each line in the data file
        for line in df:
            # Check if any value from the values file is in the line
            if any(value in line for value in values_to_filter):
                # Write the line to the output file if it contains any value
                of.write(line)

print(f"Filtered lines have been saved to {output_file}")
