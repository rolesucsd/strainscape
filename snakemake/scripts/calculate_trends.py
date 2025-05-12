#!/usr/bin/env python3
"""
Script to calculate mutation trends over time.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from scipy import stats

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def calculate_trend(
    mutations: pd.DataFrame,
    metadata: pd.DataFrame,
    p_threshold: float
) -> pd.DataFrame:
    """
    Calculate mutation trends over time.
    
    Args:
        mutations: DataFrame containing filtered mutation data
        metadata: DataFrame containing sample metadata
        p_threshold: P-value threshold for significance
        
    Returns:
        DataFrame containing trend analysis results
    """
    # Merge mutations with metadata
    merged = pd.merge(
        mutations,
        metadata,
        on='sample_id',
        how='inner'
    )
    
    # Group by mutation and calculate trend
    trends = []
    for mutation, group in merged.groupby('mutation_id'):
        # Sort by timepoint
        group = group.sort_values('timepoint')
        
        # Calculate linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            group['timepoint'],
            group['frequency']
        )
        
        # Store results
        trends.append({
            'mutation_id': mutation,
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_value ** 2,
            'p_value': p_value,
            'std_err': std_err,
            'significant': p_value < p_threshold
        })
    
    return pd.DataFrame(trends)

def main():
    """Main function."""
    # Get input and output files from Snakemake
    mutations_file = snakemake.input.mutations
    metadata_file = snakemake.input.metadata
    output_file = snakemake.output.trends
    
    # Get parameters
    p_threshold = snakemake.params.p_threshold
    
    try:
        # Load data
        logger.info("Loading mutation and metadata")
        mutations = pd.read_csv(mutations_file)
        metadata = pd.read_csv(metadata_file)
        
        # Calculate trends
        logger.info("Calculating mutation trends")
        trends = calculate_trend(
            mutations,
            metadata,
            p_threshold=p_threshold
        )
        
        # Save results
        logger.info(f"Saving trend analysis to {output_file}")
        trends.to_csv(output_file, index=False)
        
        logger.info("Trend analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error in trend analysis: {str(e)}")
        raise

if __name__ == '__main__':
    main() 