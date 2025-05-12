#!/usr/bin/env python3
"""
Script to prepare mutation trajectories over time.
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

def prepare_trajectory(
    mutations: pd.DataFrame,
    metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    Prepare mutation trajectories over time.
    
    Args:
        mutations: DataFrame containing filtered mutation data
        metadata: DataFrame containing sample metadata
        
    Returns:
        DataFrame containing trajectory analysis results
    """
    # Merge mutations with metadata
    merged = pd.merge(
        mutations,
        metadata,
        on='sample_id',
        how='inner'
    )
    
    # Sort by timepoint
    merged = merged.sort_values(['mutation_id', 'timepoint'])
    
    # Initialize results
    trajectories = []
    
    # Group by mutation
    for mutation, group in merged.groupby('mutation_id'):
        # Calculate trajectory statistics
        timepoints = group['timepoint'].values
        frequencies = group['frequency'].values
        
        # Calculate slope and intercept
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            timepoints,
            frequencies
        )
        
        # Calculate mean and std of frequency
        mean_freq = np.mean(frequencies)
        std_freq = np.std(frequencies)
        
        # Calculate frequency change
        freq_change = frequencies[-1] - frequencies[0]
        
        # Store results
        trajectories.append({
            'mutation_id': mutation,
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_value ** 2,
            'p_value': p_value,
            'std_err': std_err,
            'mean_frequency': mean_freq,
            'std_frequency': std_freq,
            'frequency_change': freq_change,
            'trajectory_type': 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'
        })
    
    return pd.DataFrame(trajectories)

def main():
    """Main function."""
    # Get input and output files from Snakemake
    mutations_file = snakemake.input.mutations
    metadata_file = snakemake.input.metadata
    output_file = snakemake.output.trajectories
    
    try:
        # Load data
        logger.info("Loading mutation data and metadata")
        mutations = pd.read_csv(mutations_file)
        metadata = pd.read_csv(metadata_file)
        
        # Prepare trajectories
        logger.info("Preparing mutation trajectories")
        trajectories = prepare_trajectory(mutations, metadata)
        
        # Save results
        logger.info(f"Saving trajectory analysis to {output_file}")
        trajectories.to_csv(output_file, index=False)
        
        logger.info("Trajectory analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error in trajectory analysis: {str(e)}")
        raise

if __name__ == '__main__':
    main() 