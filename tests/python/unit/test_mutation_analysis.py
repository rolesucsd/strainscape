#!/usr/bin/env python3
"""
Unit tests for mutation analysis functions.
"""

import pytest
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from pathlib import Path
from strainscape.analyze_mutation_types import get_mutation_type, analyze_mutation_types

@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    return pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2'],
        'Position': [3, 3, 300],
        'ref_base': ['A', 'A', 'T'],
        'new_base': ['G', 'C', 'C'],
        'coding': ['genic', 'genic', 'intergenic'],
        'Matched_Start': [1, 1, None],
        'Matched_Stop': [3, 3, None],
        'Matched_Strand': ['+', '+', None]
    })

@pytest.fixture
def sample_gene_mapping():
    """Create sample gene mapping data for testing."""
    return pd.DataFrame({
        'mutation_id': ['mut1', 'mut2', 'mut3'],
        'Matched_Gene': ['gene1', 'gene2', None],
        'Matched_Start': [90, 190, None],
        'Matched_Stop': [110, 210, None],
        'Matched_Strand': ['+', '-', None],
        'coding': ['genic', 'genic', 'intergenic']
    })

@pytest.fixture
def sample_sequences():
    """Create sample sequences for testing."""
    return {
        'chr1': 'GAA',  # GAA (E), GAG (E, silent), GAC (D, missense)
        'chr2': 'GCTAGCTAGCTAGCTAGCTA'
    }

def test_get_mutation_type():
    # Create test data
    sequences = {'chr1': 'ATGCTGCTG'}
    row = pd.Series({
        'Chromosome': 'chr1',
        'Position': 4,
        'ref_base': 'C',
        'new_base': 'T',
        'gene_type': 'genic',
        'Start': 1,
        'Stop': 9,
        'Strand': '+'
    })
    
    # Test the function
    mut_type, coding = get_mutation_type(row, sequences)
    assert mut_type in ['Silent', 'Missense', 'Nonsense']
    assert coding in ['Coding', 'Non-Coding', 'Error']

def test_analyze_mutations():
    # Create test data
    sequences = {'chr1': 'ATGCTGCTG'}
    mutations = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Position': [4, 7],
        'ref_base': ['C', 'G'],
        'new_base': ['T', 'A'],
        'gene_type': ['genic', 'genic'],
        'Start': [1, 1],
        'Stop': [9, 9],
        'Strand': ['+', '+']
    })
    
    # Test the function
    result = analyze_mutation_types(mutations, sequences)
    assert 'Mutation_Type' in result.columns
    assert 'Coding_Status' in result.columns
    assert len(result) == len(mutations)

def test_analyze_mutations_empty_input():
    """Test mutation analysis with empty input."""
    # Create empty DataFrame
    empty_df = pd.DataFrame(columns=['Chromosome', 'Position', 'ref_base', 'new_base', 
                                   'coding', 'Matched_Start', 'Matched_Stop', 'Matched_Strand'])
    
    # Run analysis
    results = analyze_mutation_types(empty_df, {})
    
    # Check results
    assert len(results) == 0
    assert 'Mutation_Type' in results.columns 