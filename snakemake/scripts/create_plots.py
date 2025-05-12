#!/usr/bin/env python3
"""
Script to create visualization plots for analysis results.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_mutation_frequency_plot(
    mutations: pd.DataFrame,
    metadata: pd.DataFrame,
    pdf: PdfPages
) -> None:
    """
    Create plot of mutation frequencies over time.
    
    Args:
        mutations: DataFrame containing filtered mutation data
        metadata: DataFrame containing sample metadata
        pdf: PDF file to save plot
    """
    # Merge data
    merged = pd.merge(
        mutations,
        metadata,
        on='sample_id',
        how='inner'
    )
    
    # Create plot
    plt.figure(figsize=(12, 8))
    sns.lineplot(
        data=merged,
        x='timepoint',
        y='frequency',
        hue='mutation_id',
        marker='o'
    )
    plt.title('Mutation Frequencies Over Time')
    plt.xlabel('Timepoint')
    plt.ylabel('Frequency')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    pdf.savefig()
    plt.close()

def create_trend_plot(
    trends: pd.DataFrame,
    pdf: PdfPages
) -> None:
    """
    Create plot of mutation trends.
    
    Args:
        trends: DataFrame containing trend analysis results
        pdf: PDF file to save plot
    """
    # Create plot
    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=trends,
        x='slope',
        hue='significant',
        multiple='stack'
    )
    plt.title('Distribution of Mutation Trends')
    plt.xlabel('Slope')
    plt.ylabel('Count')
    plt.tight_layout()
    
    # Save plot
    pdf.savefig()
    plt.close()

def create_gene_analysis_plot(
    analysis: pd.DataFrame,
    pdf: PdfPages
) -> None:
    """
    Create plot of gene analysis results.
    
    Args:
        analysis: DataFrame containing mutation analysis results
        pdf: PDF file to save plot
    """
    # Create plot
    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=analysis,
        x='gene_id',
        y='total_mutations'
    )
    plt.title('Total Mutations per Gene')
    plt.xlabel('Gene')
    plt.ylabel('Number of Mutations')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save plot
    pdf.savefig()
    plt.close()

def create_trajectory_plot(
    trajectories: pd.DataFrame,
    pdf: PdfPages
) -> None:
    """
    Create plot of mutation trajectories.
    
    Args:
        trajectories: DataFrame containing trajectory analysis results
        pdf: PDF file to save plot
    """
    # Create plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=trajectories,
        x='trajectory_type',
        y='frequency_change'
    )
    plt.title('Frequency Change by Trajectory Type')
    plt.xlabel('Trajectory Type')
    plt.ylabel('Frequency Change')
    plt.tight_layout()
    
    # Save plot
    pdf.savefig()
    plt.close()

def main():
    """Main function."""
    # Get input and output files from Snakemake
    mutations_file = snakemake.input.mutations
    trends_file = snakemake.input.trends
    mapping_file = snakemake.input.mapping
    analysis_file = snakemake.input.analysis
    trajectories_file = snakemake.input.trajectories
    output_file = snakemake.output.plots
    
    try:
        # Load data
        logger.info("Loading analysis results")
        mutations = pd.read_csv(mutations_file)
        trends = pd.read_csv(trends_file)
        mapping = pd.read_csv(mapping_file)
        analysis = pd.read_csv(analysis_file)
        trajectories = pd.read_csv(trajectories_file)
        
        # Create plots
        logger.info("Creating plots")
        with PdfPages(output_file) as pdf:
            # Create mutation frequency plot
            create_mutation_frequency_plot(mutations, mapping, pdf)
            
            # Create trend plot
            create_trend_plot(trends, pdf)
            
            # Create gene analysis plot
            create_gene_analysis_plot(analysis, pdf)
            
            # Create trajectory plot
            create_trajectory_plot(trajectories, pdf)
        
        logger.info("Plot creation completed successfully")
        
    except Exception as e:
        logger.error(f"Error in plot creation: {str(e)}")
        raise

if __name__ == '__main__':
    main() 