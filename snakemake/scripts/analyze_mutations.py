#!/usr/bin/env python3
"""
Script to analyze mutation types and patterns.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from collections import Counter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_mutation_type(ref: str, alt: str) -> str:
    """
    Determine mutation type from reference and alternate alleles.
    
    Args:
        ref: Reference allele
        alt: Alternate allele
        
    Returns:
        Mutation type (e.g., 'SNP', 'INDEL')
    """
    if len(ref) == len(alt) == 1:
        return 'SNP'
    elif len(ref) != len(alt):
        return 'INDEL'
    else:
        return 'COMPLEX'

def analyze_mutations(
    mutations: pd.DataFrame,
    mapping: pd.DataFrame
) -> pd.DataFrame:
    """
    Analyze mutation types and patterns.
    
    Args:
        mutations: DataFrame containing filtered mutation data
        mapping: DataFrame containing gene mapping results
        
    Returns:
        DataFrame containing mutation analysis results
    """
    # Merge mutations with gene mapping
    merged = pd.merge(
        mutations,
        mapping,
        on='mutation_id',
        how='inner'
    )
    
    # Initialize results
    analysis = []
    
    # Group by gene
    for gene, group in merged.groupby('gene_id'):
        # Count mutation types
        mutation_types = Counter(
            get_mutation_type(row['ref'], row['alt'])
            for _, row in group.iterrows()
        )
        
        # Calculate statistics
        total_mutations = len(group)
        snp_count = mutation_types.get('SNP', 0)
        indel_count = mutation_types.get('INDEL', 0)
        complex_count = mutation_types.get('COMPLEX', 0)
        
        # Store results
        analysis.append({
            'gene_id': gene,
            'total_mutations': total_mutations,
            'snp_count': snp_count,
            'indel_count': indel_count,
            'complex_count': complex_count,
            'snp_fraction': snp_count / total_mutations if total_mutations > 0 else 0,
            'indel_fraction': indel_count / total_mutations if total_mutations > 0 else 0,
            'complex_fraction': complex_count / total_mutations if total_mutations > 0 else 0
        })
    
    return pd.DataFrame(analysis)

def main():
    """Main function."""
    # Get input and output files from Snakemake
    mutations_file = snakemake.input.mutations
    mapping_file = snakemake.input.mapping
    output_file = snakemake.output.analysis
    
    try:
        # Load data
        logger.info("Loading mutation data and gene mapping")
        mutations = pd.read_csv(mutations_file)
        mapping = pd.read_csv(mapping_file)
        
        # Analyze mutations
        logger.info("Analyzing mutations")
        analysis = analyze_mutations(mutations, mapping)
        
        # Save results
        logger.info(f"Saving mutation analysis to {output_file}")
        analysis.to_csv(output_file, index=False)
        
        logger.info("Mutation analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error in mutation analysis: {str(e)}")
        raise

if __name__ == '__main__':
    main() 