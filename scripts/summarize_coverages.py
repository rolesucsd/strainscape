import os
import pandas as pd
from glob import glob

# Directory containing the coverage results
coverage_dir = 'coverage_results'

# Use glob to find all coverage files
coverage_files = glob(os.path.join(coverage_dir, '*_coverage.txt'))

# Initialize an empty list to store the summarized results
summary_data = []

# Process each file
for file_path in coverage_files:
    sample_name = os.path.basename(file_path).replace('_coverage.txt', '')
    print(sample_name)
    try:
        # Load the coverage file
        coverage_df = pd.read_csv(file_path, sep='\t', header=None, names=['contig', 'position', 'coverage'])
        # Calculate average coverage
        average_coverage = coverage_df['coverage'].mean() if not coverage_df.empty else 0
        # Calculate the percentage of positions covered
        total_positions = len(coverage_df)
        positions_covered = (coverage_df['coverage'] > 0).sum()
        coverage_percent = (positions_covered / total_positions * 100) if total_positions else 0
        summary_data.append({
            'Sample': sample_name,
            'Average_Coverage': average_coverage,
            'Coverage_Percent': coverage_percent
        })
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

# Create a DataFrame from the list of dictionaries
summary_df = pd.DataFrame(summary_data)

# Save the summary to a file
summary_df.to_csv('coverage_summary.csv', index=False)

print(f"Processed {len(coverage_files)} files.")
