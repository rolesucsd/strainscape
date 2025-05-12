"""Module for analyzing and summarizing coverage data."""

import os
import pandas as pd
from glob import glob

def process_coverage_file(file_path):
    """Process a single coverage file and return summary statistics.
    
    Args:
        file_path (str): Path to the coverage file
        
    Returns:
        dict: Dictionary containing sample name and coverage statistics
    """
    sample_name = os.path.basename(file_path).replace('_coverage.txt', '')
    try:
        coverage_df = pd.read_csv(file_path, sep='\t', header=None, 
                                names=['contig', 'position', 'coverage'])
        average_coverage = coverage_df['coverage'].mean() if not coverage_df.empty else 0
        total_positions = len(coverage_df)
        positions_covered = (coverage_df['coverage'] > 0).sum()
        coverage_percent = (positions_covered / total_positions * 100) if total_positions else 0
        
        return {
            'Sample': sample_name,
            'Average_Coverage': average_coverage,
            'Coverage_Percent': coverage_percent
        }
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def summarize_coverage(coverage_dir, output_file='coverage_summary.csv'):
    """Summarize coverage data from multiple files.
    
    Args:
        coverage_dir (str): Directory containing coverage files
        output_file (str): Path to save the summary CSV file
        
    Returns:
        pd.DataFrame: DataFrame containing coverage summaries
    """
    coverage_files = glob(os.path.join(coverage_dir, '*_coverage.txt'))
    summary_data = []
    
    for file_path in coverage_files:
        result = process_coverage_file(file_path)
        if result:
            summary_data.append(result)
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False)
    print(f"Processed {len(coverage_files)} files.")
    return summary_df

def main():
    """Command-line interface for summarizing coverage data."""
    import argparse
    parser = argparse.ArgumentParser(description="Summarize coverage data from multiple files.")
    parser.add_argument("--coverage_dir", required=True, help="Directory containing coverage files")
    parser.add_argument("--output_file", default="coverage_summary.csv", 
                       help="Output file for coverage summary")
    args = parser.parse_args()
    
    summarize_coverage(args.coverage_dir, args.output_file)

if __name__ == "__main__":
    main() 