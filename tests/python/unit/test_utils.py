#!/usr/bin/env python3
"""
Unit tests for utility functions.
"""

import pytest
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import tempfile
import os

# Add the directory containing the Snakemake scripts to the Python path
sys.path.append(str(Path(__file__).parent.parent.parent / "snakemake" / "scripts"))

# Import the functions to test
from utils import (
    load_data,
    save_data,
    filter_data,
    merge_data,
    calculate_statistics
)

@pytest.fixture
def sample_data():
    """Create sample data for testing."""
    return pd.DataFrame({
        'id': ['1', '2', '3', '4', '5'],
        'value': [1.0, 2.0, 3.0, 4.0, 5.0],
        'category': ['A', 'B', 'A', 'B', 'A']
    })

@pytest.fixture
def sample_data2():
    """Create another sample dataset for testing."""
    return pd.DataFrame({
        'id': ['1', '2', '3', '4', '5'],
        'score': [10, 20, 30, 40, 50],
        'category': ['A', 'B', 'A', 'B', 'A']
    })

def test_load_data():
    """Test loading data from CSV file."""
    # Create temporary CSV file
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as temp:
        temp.write(b'id,value,category\n1,1.0,A\n2,2.0,B\n3,3.0,A\n')
        temp_path = temp.name
    
    try:
        # Load data
        data = load_data(temp_path)
        
        # Check structure
        assert isinstance(data, pd.DataFrame)
        assert len(data) == 3
        assert all(col in data.columns for col in ['id', 'value', 'category'])
        
        # Check content
        assert data.iloc[0]['id'] == '1'
        assert data.iloc[0]['value'] == 1.0
        assert data.iloc[0]['category'] == 'A'
    finally:
        os.unlink(temp_path)

def test_save_data(sample_data):
    """Test saving data to CSV file."""
    # Create temporary file path
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as temp:
        temp_path = temp.name
    
    try:
        # Save data
        save_data(sample_data, temp_path)
        
        # Load saved data
        loaded_data = pd.read_csv(temp_path)
        
        # Check content
        assert len(loaded_data) == len(sample_data)
        assert all(col in loaded_data.columns for col in sample_data.columns)
        assert loaded_data.equals(sample_data)
    finally:
        os.unlink(temp_path)

def test_filter_data(sample_data):
    """Test filtering data based on conditions."""
    # Filter by value
    filtered = filter_data(sample_data, 'value > 2.0')
    assert len(filtered) == 3
    assert all(filtered['value'] > 2.0)
    
    # Filter by category
    filtered = filter_data(sample_data, "category == 'A'")
    assert len(filtered) == 3
    assert all(filtered['category'] == 'A')
    
    # Filter with multiple conditions
    filtered = filter_data(sample_data, "value > 2.0 and category == 'A'")
    assert len(filtered) == 2
    assert all(filtered['value'] > 2.0)
    assert all(filtered['category'] == 'A')

def test_merge_data(sample_data, sample_data2):
    """Test merging two datasets."""
    # Merge on id
    merged = merge_data(sample_data, sample_data2, on='id')
    
    # Check structure
    assert isinstance(merged, pd.DataFrame)
    assert len(merged) == 5
    assert all(col in merged.columns for col in ['id', 'value', 'score', 'category'])
    
    # Check content
    assert merged.iloc[0]['id'] == '1'
    assert merged.iloc[0]['value'] == 1.0
    assert merged.iloc[0]['score'] == 10
    
    # Test with different merge types
    merged_left = merge_data(sample_data, sample_data2, on='id', how='left')
    assert len(merged_left) == 5
    
    merged_right = merge_data(sample_data, sample_data2, on='id', how='right')
    assert len(merged_right) == 5
    
    merged_outer = merge_data(sample_data, sample_data2, on='id', how='outer')
    assert len(merged_outer) == 5

def test_calculate_statistics(sample_data):
    """Test calculating basic statistics."""
    # Calculate statistics for numeric column
    stats = calculate_statistics(sample_data, 'value')
    
    # Check structure
    assert isinstance(stats, dict)
    assert 'mean' in stats
    assert 'std' in stats
    assert 'min' in stats
    assert 'max' in stats
    
    # Check values
    assert stats['mean'] == 3.0
    assert stats['std'] == pytest.approx(1.5811, rel=1e-4)
    assert stats['min'] == 1.0
    assert stats['max'] == 5.0
    
    # Test with groupby
    group_stats = calculate_statistics(sample_data, 'value', groupby='category')
    assert isinstance(group_stats, pd.DataFrame)
    assert len(group_stats) == 2
    assert all(col in group_stats.columns for col in ['category', 'mean', 'std', 'min', 'max'])

def test_load_data_invalid_file():
    """Test loading data from invalid file."""
    with pytest.raises(FileNotFoundError):
        load_data('nonexistent_file.csv')

def test_save_data_invalid_path(sample_data):
    """Test saving data to invalid path."""
    with pytest.raises(PermissionError):
        save_data(sample_data, '/nonexistent/path/data.csv')

def test_filter_data_invalid_condition(sample_data):
    """Test filtering data with invalid condition."""
    with pytest.raises(Exception):
        filter_data(sample_data, 'invalid_condition')

def test_merge_data_invalid_column(sample_data, sample_data2):
    """Test merging data with invalid column."""
    with pytest.raises(KeyError):
        merge_data(sample_data, sample_data2, on='nonexistent_column')

def test_calculate_statistics_invalid_column(sample_data):
    """Test calculating statistics for invalid column."""
    with pytest.raises(KeyError):
        calculate_statistics(sample_data, 'nonexistent_column') 