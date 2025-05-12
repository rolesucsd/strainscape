#!/usr/bin/env python3
"""
Script to map mutations to genes using phylogenetic tree.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_tree(tree_file: str) -> Phylo.BaseTree.Tree:
    """
    Read phylogenetic tree from file.
    
    Args:
        tree_file: Path to tree file
        
    Returns:
        Phylogenetic tree object
    """
    return Phylo.read(tree_file, 'newick')

def map_mutations_to_genes(
    mutations: pd.DataFrame,
    tree: Phylo.BaseTree.Tree
) -> pd.DataFrame:
    """
    Map mutations to genes using phylogenetic tree.
    
    Args:
        mutations: DataFrame containing filtered mutation data
        tree: Phylogenetic tree object
        
    Returns:
        DataFrame containing gene mapping results
    """
    # Initialize results
    mapping = []
    
    # Get all terminal nodes (genes)
    genes = [node.name for node in tree.get_terminals()]
    
    # For each mutation
    for mutation, group in mutations.groupby('mutation_id'):
        # Get mutation position
        position = group['position'].iloc[0]
        
        # Find closest gene
        closest_gene = None
        min_distance = float('inf')
        
        for gene in genes:
            # Get gene position from tree
            gene_node = tree.find_any(gene)
            gene_position = gene_node.branch_length
            
            # Calculate distance
            distance = abs(position - gene_position)
            
            if distance < min_distance:
                min_distance = distance
                closest_gene = gene
        
        # Store mapping
        mapping.append({
            'mutation_id': mutation,
            'gene_id': closest_gene,
            'distance': min_distance
        })
    
    return pd.DataFrame(mapping)

def main():
    """Main function."""
    # Get input and output files from Snakemake
    mutations_file = snakemake.input.mutations
    tree_file = snakemake.input.tree
    output_file = snakemake.output.mapping
    
    try:
        # Load data
        logger.info("Loading mutation data and tree")
        mutations = pd.read_csv(mutations_file)
        tree = read_tree(tree_file)
        
        # Map mutations to genes
        logger.info("Mapping mutations to genes")
        mapping = map_mutations_to_genes(mutations, tree)
        
        # Save results
        logger.info(f"Saving gene mapping to {output_file}")
        mapping.to_csv(output_file, index=False)
        
        logger.info("Gene mapping completed successfully")
        
    except Exception as e:
        logger.error(f"Error in gene mapping: {str(e)}")
        raise

if __name__ == '__main__':
    main() 