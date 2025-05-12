#!/usr/bin/env python3
"""
Script to filter mutations based on coverage and frequency thresholds.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def filter_mutations(
    mutation_data: pd.DataFrame,
    min_coverage: int,
    min_freq: float
) -> pd.DataFrame:
    """
    Filter mutations based on coverage and frequency thresholds.
    
    Args:
        mutation_data: DataFrame containing mutation data
        min_coverage: Minimum coverage threshold
        min_freq: Minimum frequency threshold
        
    Returns:
        Filtered DataFrame
    """
    # Filter by coverage
    filtered = mutation_data[mutation_data['coverage'] >= min_coverage]
    
    # Filter by frequency
    filtered = filtered[filtered['frequency'] >= min_freq]
    
    return filtered

def main():
    """Main function."""
    # Get input and output files from Snakemake
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    
    # Get parameters
    min_coverage = snakemake.params.min_coverage
    min_freq = snakemake.params.min_freq
    
    try:
        # Load data
        logger.info(f"Loading mutation data from {input_file}")
        mutation_data = pd.read_csv(input_file)
        
        # Filter mutations
        logger.info("Filtering mutations")
        filtered = filter_mutations(
            mutation_data,
            min_coverage=min_coverage,
            min_freq=min_freq
        )
        
        # Save results
        logger.info(f"Saving filtered mutations to {output_file}")
        filtered.to_csv(output_file, index=False)
        
        logger.info("Filtering completed successfully")
        
    except Exception as e:
        logger.error(f"Error in filtering: {str(e)}")
        raise

if __name__ == '__main__':
    main() 