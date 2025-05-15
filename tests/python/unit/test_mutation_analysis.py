#!/usr/bin/env python3
"""
Unit tests for mutation analysis functions.
"""

import pytest
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import sys
from pathlib import Path

# Add the directory containing the Snakemake scripts to the Python path
sys.path.append(str(Path(__file__).parent.parent.parent / "snakemake" / "scripts"))

# Import the functions to test
from analyze_mutation_types import get_mutation_type, load_sequences

@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    return pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2'],
        'Position': [100, 200, 300],
        'ref_base': ['A', 'G', 'T'],
        'new_base': ['G', 'A', 'C'],
        'coding': ['genic', 'genic', 'intergenic'],
        'Matched_Start': [90, 190, None],
        'Matched_Stop': [110, 210, None],
        'Matched_Strand': ['+', '-', None]
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
        'chr1': 'ATGCATGCATGCATGCATGC',  # 20 bases
        'chr2': 'GCTAGCTAGCTAGCTAGCTA'   # 20 bases
    }

def test_filter_mutations():
    """Test mutation filtering function."""
    from filter_mutations import filter_mutations
    
    # Create sample data
    mutations = pd.DataFrame({
        'mutation_id': ['mut1', 'mut2', 'mut3'],
        'coverage': [100, 50, 200],
        'frequency': [0.8, 0.3, 0.9]
    })
    
    # Test filtering
    filtered = filter_mutations(mutations, min_coverage=60, min_freq=0.4)
    
    # Check results
    assert len(filtered) == 2
    assert 'mut1' in filtered['mutation_id'].values
    assert 'mut3' in filtered['mutation_id'].values
    assert 'mut2' not in filtered['mutation_id'].values

def test_get_mutation_type(sample_mutation_data, sample_sequences):
    """Test mutation type determination."""
    # Test silent mutation
    row = sample_mutation_data.iloc[0]
    result = get_mutation_type(row, sample_sequences)
    assert result['Mutation_Type'] in ['Silent', 'Missense', 'Nonsense']
    
    # Test intergenic mutation
    row = sample_mutation_data.iloc[2]
    result = get_mutation_type(row, sample_sequences)
    assert result['Mutation_Type'] is None

def test_analyze_mutations(sample_mutation_data, sample_sequences):
    """Test mutation analysis function."""
    from analyze_mutation_types import analyze_mutation_types
    
    # Run analysis
    results = analyze_mutation_types(sample_mutation_data, sample_sequences)
    
    # Check results structure
    assert 'Mutation_Type' in results.columns
    assert len(results) == len(sample_mutation_data)
    
    # Check mutation type counts
    type_counts = results['Mutation_Type'].value_counts()
    assert type_counts.get('Silent', 0) + type_counts.get('Missense', 0) + type_counts.get('Nonsense', 0) == 2
    assert type_counts.get(None, 0) == 1

def test_analyze_mutations_empty_input():
    """Test mutation analysis with empty input."""
    from analyze_mutation_types import analyze_mutation_types
    
    # Create empty DataFrame
    empty_df = pd.DataFrame(columns=['Chromosome', 'Position', 'ref_base', 'new_base', 
                                   'coding', 'Matched_Start', 'Matched_Stop', 'Matched_Strand'])
    
    # Run analysis
    results = analyze_mutation_types(empty_df, {})
    
    # Check results
    assert len(results) == 0
    assert 'Mutation_Type' in results.columns 