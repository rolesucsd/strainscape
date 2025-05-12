import os
import pandas as pd

# Directory containing the coverage results
coverage_dir = 'coverage_results'

# Assuming amrfinder_df is loaded here with the necessary columns ['Gene symbol', 'Contig id', 'Start', 'Stop']
# Load or define amrfinder_df here
amrfinder_df = pd.read_csv('amrfinder_output.txt', sep="\t")

# Initialize an empty list to store the summarized results
summary_data = []

# Loop through each coverage file
for coverage_file in os.listdir(coverage_dir):
    if coverage_file.endswith('_coverage.txt'):
        sample_name = coverage_file.replace('_coverage.txt', '')
        
        # Load the coverage file
        file_path = os.path.join(coverage_dir, coverage_file)
        try:
            coverage_df = pd.read_csv(file_path, sep='\t', header=None, names=['contig', 'position', 'coverage'])
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
        
        # Get unique genes to be analyzed
        unique_genes = amrfinder_df['Gene symbol'].unique()
        
        # Summarize the coverage by gene
        for gene_name in unique_genes:
            gene_info = amrfinder_df[amrfinder_df['Gene symbol'] == gene_name].iloc[0]
            contig_id = gene_info['Gene symbol']
            start = 1
            stop = gene_info['Stop'] - gene_info['Start']
            
            gene_coverage = coverage_df[(coverage_df['contig'] == contig_id) & 
                                        (coverage_df['position'] >= start) & 
                                        (coverage_df['position'] <= stop)]
            
            print(gene_coverage)
            
            average_coverage = gene_coverage['coverage'].mean() if not gene_coverage.empty else 0
            
            summary_data.append({
                'Sample': sample_name,
                'Gene': gene_name,
                'Average_Coverage': average_coverage
            })

# Convert the list of dicts to a DataFrame
summary_df = pd.DataFrame(summary_data)

# Save the summary to a file
summary_df.to_csv('coverage_summary.csv', index=False)
