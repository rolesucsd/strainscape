"""Tests for the trees module."""

import os
import tempfile
import pandas as pd
import pytest
from strainscape.trees import read_annotation, build_and_save_tree, compute_and_save_trees

@pytest.fixture
def sample_genes():
    """Create a sample genes DataFrame for testing."""
    return pd.DataFrame({
        'Sequence Id': ['chr1', 'chr1', 'chr2'],
        'Start': [100, 200, 300],
        'Stop': [150, 250, 350],
        'Locus Tag': ['gene1', 'gene2', 'gene3']
    })

def test_read_annotation(sample_genes, tmp_path):
    """Test reading gene annotations."""
    # Save sample data to a temporary file
    genes_file = tmp_path / "test_genes.txt"
    sample_genes.to_csv(genes_file, sep='\t', index=False)
    
    # Read the file
    result = read_annotation(genes_file)
    
    # Check results
    assert 'Chromosome' in result.columns
    assert 'Start' in result.columns
    assert 'Stop' in result.columns
    assert len(result) == 3
    assert result['Start'].dtype == 'int64'
    assert result['Stop'].dtype == 'int64'

def test_build_and_save_tree(sample_genes, tmp_path):
    """Test building and saving an interval tree."""
    # Group by chromosome
    chrom = 'chr1'
    chrom_df = sample_genes[sample_genes['Sequence Id'] == chrom]
    
    # Build and save tree
    result_chrom, result_tree = build_and_save_tree(chrom, chrom_df, tmp_path)
    
    # Check results
    assert result_chrom == chrom
    assert len(result_tree) == 2  # Two genes in chr1
    
    # Check if file was created
    tree_file = tmp_path / f"chromosome_tree_{chrom}.pkl"
    assert tree_file.exists()

def test_compute_and_save_trees(sample_genes, tmp_path):
    """Test computing and saving multiple trees."""
    # Save sample data
    genes_file = tmp_path / "test_genes.txt"
    sample_genes.to_csv(genes_file, sep='\t', index=False)
    
    # Compute trees
    results = compute_and_save_trees(genes_file, tmp_path)
    
    # Check results
    assert 'chr1' in results
    assert 'chr2' in results
    assert len(results['chr1']) == 2
    assert len(results['chr2']) == 1 