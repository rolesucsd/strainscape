#!/usr/bin/env python3
"""
Unit tests for trajectory analysis functions.
"""

import pytest
import pandas as pd
import numpy as np
from scipy import stats
import sys
from pathlib import Path

# Add the directory containing the Snakemake scripts to the Python path
sys.path.append(str(Path(__file__).parent.parent.parent / "snakemake" / "scripts"))

# Import the functions to test
from prepare_trajectories import filter_static_mutations, prepare_trajectories

@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    return pd.DataFrame({
        'mutation_id': ['mut1', 'mut1', 'mut1', 'mut2', 'mut2', 'mut2', 'mut3', 'mut3', 'mut3'],
        'timepoint': [0, 1, 2, 0, 1, 2, 0, 1, 2],
        'frequency': [0.1, 0.2, 0.3, 0.5, 0.5, 0.5, 0.8, 0.6, 0.4]  # mut1: increasing, mut2: static, mut3: decreasing
    })

@pytest.fixture
def sample_metadata():
    """Create sample metadata for testing."""
    return pd.DataFrame({
        'timepoint': [0, 1, 2],
        'sample_id': ['sample1', 'sample2', 'sample3']
    })

def test_filter_static_mutations(sample_mutation_data):
    """Test filtering of static mutations."""
    # Test with different thresholds
    filtered = filter_static_mutations(sample_mutation_data, min_freq_change=0.1)
    assert len(filtered['mutation_id'].unique()) == 2
    assert 'mut1' in filtered['mutation_id'].values
    assert 'mut3' in filtered['mutation_id'].values
    assert 'mut2' not in filtered['mutation_id'].values
    
    # Test with higher threshold
    filtered = filter_static_mutations(sample_mutation_data, min_freq_change=0.3)
    assert len(filtered['mutation_id'].unique()) == 1
    assert 'mut3' in filtered['mutation_id'].values

def test_prepare_trajectories(sample_mutation_data, sample_metadata):
    """Test trajectory preparation."""
    # Prepare trajectories
    trajectories = prepare_trajectories(sample_mutation_data, sample_metadata, p_threshold=0.05)
    
    # Check results structure
    assert isinstance(trajectories, pd.DataFrame)
    assert 'mutation_id' in trajectories.columns
    assert 'slope' in trajectories.columns
    assert 'intercept' in trajectories.columns
    assert 'r_squared' in trajectories.columns
    assert 'p_value' in trajectories.columns
    assert 'mean_frequency' in trajectories.columns
    assert 'std_frequency' in trajectories.columns
    assert 'frequency_change' in trajectories.columns
    assert 'trajectory_type' in trajectories.columns
    
    # Check specific results
    mut1_traj = trajectories[trajectories['mutation_id'] == 'mut1'].iloc[0]
    assert mut1_traj['slope'] > 0
    assert mut1_traj['trajectory_type'] == 'increasing'
    assert mut1_traj['frequency_change'] == 0.2  # 0.3 - 0.1
    
    mut2_traj = trajectories[trajectories['mutation_id'] == 'mut2'].iloc[0]
    assert mut2_traj['slope'] == 0
    assert mut2_traj['trajectory_type'] == 'static'
    assert mut2_traj['frequency_change'] == 0.0  # 0.5 - 0.5
    
    mut3_traj = trajectories[trajectories['mutation_id'] == 'mut3'].iloc[0]
    assert mut3_traj['slope'] < 0
    assert mut3_traj['trajectory_type'] == 'decreasing'
    assert mut3_traj['frequency_change'] == -0.4  # 0.4 - 0.8

def test_prepare_trajectories_small_sample():
    """Test trajectory preparation with small sample size."""
    # Create small sample data
    data = pd.DataFrame({
        'mutation_id': ['mut1', 'mut1'],
        'timepoint': [0, 1],
        'frequency': [0.1, 0.9]
    })
    metadata = pd.DataFrame({
        'timepoint': [0, 1],
        'sample_id': ['sample1', 'sample2']
    })
    
    # Prepare trajectories
    trajectories = prepare_trajectories(data, metadata, p_threshold=0.05)
    
    # Check results
    assert len(trajectories) == 1
    assert trajectories.iloc[0]['mutation_id'] == 'mut1'
    assert trajectories.iloc[0]['frequency_change'] == 0.8
    assert trajectories.iloc[0]['trajectory_type'] == 'increasing'
    assert 'p_value' in trajectories.columns
    assert 'significant' in trajectories.columns

def test_prepare_trajectories_empty_input():
    """Test trajectory preparation with empty input."""
    empty_mutations = pd.DataFrame(columns=['mutation_id', 'timepoint', 'frequency'])
    empty_metadata = pd.DataFrame(columns=['timepoint', 'sample_id'])
    
    trajectories = prepare_trajectories(empty_mutations, empty_metadata, p_threshold=0.05)
    
    assert len(trajectories) == 0
    assert all(col in trajectories.columns for col in [
        'mutation_id', 'slope', 'intercept', 'r_squared',
        'p_value', 'mean_frequency', 'std_frequency',
        'frequency_change', 'trajectory_type'
    ]) 