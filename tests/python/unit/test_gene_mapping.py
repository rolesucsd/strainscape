#!/usr/bin/env python3
"""
Unit tests for gene mapping functions.
"""

import pytest
import pandas as pd
import numpy as np
import pickle
import sys
from pathlib import Path
from Bio import SeqIO
from io import StringIO

# Add the directory containing the Snakemake scripts to the Python path
sys.path.append(str(Path(__file__).parent.parent.parent / "snakemake" / "scripts"))

# Import the functions to test
from map_genes import load_gene_tree, map_mutations_to_genes

@pytest.fixture
def sample_gene_tree():
    """Create a sample gene tree for testing."""
    # Create a simple tree structure
    tree = {
        'gene1': {
            'mutations': ['mut1', 'mut2'],
            'children': {
                'gene1.1': {'mutations': ['mut3']},
                'gene1.2': {'mutations': ['mut4']}
            }
        },
        'gene2': {
            'mutations': ['mut5'],
            'children': {}
        }
    }
    return tree

@pytest.fixture
def sample_mutations():
    """Create sample mutation data for testing."""
    return pd.DataFrame({
        'mutation_id': ['mut1', 'mut2', 'mut3', 'mut4', 'mut5', 'mut6'],
        'position': [100, 200, 300, 400, 500, 600],
        'ref_base': ['A', 'C', 'G', 'T', 'A', 'C'],
        'new_base': ['G', 'T', 'A', 'C', 'G', 'T']
    })

@pytest.fixture
def sample_fasta():
    """Create a sample FASTA file content for testing."""
    return StringIO(""">gene1
ATCGATCG
>gene2
GCTAGCTA
""")

def test_load_gene_tree(sample_gene_tree, tmp_path):
    """Test loading gene tree from pickle file."""
    # Save sample tree to temporary file
    tree_file = tmp_path / "test_tree.pkl"
    with open(tree_file, 'wb') as f:
        pickle.dump(sample_gene_tree, f)
    
    # Load tree
    loaded_tree = load_gene_tree(str(tree_file))
    
    # Check structure
    assert isinstance(loaded_tree, dict)
    assert 'gene1' in loaded_tree
    assert 'gene2' in loaded_tree
    assert 'mutations' in loaded_tree['gene1']
    assert 'children' in loaded_tree['gene1']
    
    # Check content
    assert loaded_tree['gene1']['mutations'] == ['mut1', 'mut2']
    assert loaded_tree['gene2']['mutations'] == ['mut5']
    assert 'gene1.1' in loaded_tree['gene1']['children']
    assert 'gene1.2' in loaded_tree['gene1']['children']

def test_map_mutations_to_genes(sample_mutations, sample_gene_tree, sample_fasta):
    """Test mapping mutations to genes."""
    # Create temporary FASTA file
    fasta_file = StringIO()
    SeqIO.write(SeqIO.parse(sample_fasta, "fasta"), fasta_file, "fasta")
    fasta_file.seek(0)
    
    # Map mutations
    mapped_mutations = map_mutations_to_genes(sample_mutations, sample_gene_tree, fasta_file)
    
    # Check results structure
    assert isinstance(mapped_mutations, pd.DataFrame)
    assert 'mutation_id' in mapped_mutations.columns
    assert 'gene' in mapped_mutations.columns
    assert 'mutation_type' in mapped_mutations.columns
    
    # Check specific mappings
    mut1_mapping = mapped_mutations[mapped_mutations['mutation_id'] == 'mut1'].iloc[0]
    assert mut1_mapping['gene'] == 'gene1'
    assert mut1_mapping['mutation_type'] in ['silent', 'missense', 'nonsense']
    
    mut5_mapping = mapped_mutations[mapped_mutations['mutation_id'] == 'mut5'].iloc[0]
    assert mut5_mapping['gene'] == 'gene2'
    assert mut5_mapping['mutation_type'] in ['silent', 'missense', 'nonsense']
    
    # Check unmapped mutations
    mut6_mapping = mapped_mutations[mapped_mutations['mutation_id'] == 'mut6']
    assert len(mut6_mapping) == 0

def test_map_mutations_to_genes_empty_input():
    """Test mapping mutations with empty input."""
    empty_mutations = pd.DataFrame(columns=['mutation_id', 'position', 'ref_base', 'new_base'])
    empty_tree = {}
    empty_fasta = StringIO("")
    
    mapped_mutations = map_mutations_to_genes(empty_mutations, empty_tree, empty_fasta)
    
    assert len(mapped_mutations) == 0
    assert all(col in mapped_mutations.columns for col in ['mutation_id', 'gene', 'mutation_type'])

def test_map_mutations_to_genes_invalid_tree():
    """Test mapping mutations with invalid tree structure."""
    mutations = pd.DataFrame({
        'mutation_id': ['mut1'],
        'position': [100],
        'ref_base': ['A'],
        'new_base': ['G']
    })
    invalid_tree = {'gene1': 'invalid_structure'}
    fasta_file = StringIO(">gene1\nATCG")
    
    with pytest.raises(KeyError):
        map_mutations_to_genes(mutations, invalid_tree, fasta_file) 