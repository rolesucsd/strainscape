#!/usr/bin/env python3
"""
Unit tests for trend analysis functions.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os

# Add the snakemake scripts directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../snakemake/scripts'))

from calculate_trends import calculate_trend
from prepare_trajectories import prepare_trajectory

@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    return pd.DataFrame({
        'mutation_id': ['mut1', 'mut1', 'mut2', 'mut2'],
        'sample_id': ['sample1', 'sample2', 'sample1', 'sample2'],
        'frequency': [0.2, 0.8, 0.9, 0.1]
    })

@pytest.fixture
def sample_metadata():
    """Create sample metadata for testing."""
    return pd.DataFrame({
        'sample_id': ['sample1', 'sample2'],
        'timepoint': [0, 10]
    })

def test_calculate_trend(sample_mutation_data, sample_metadata):
    """Test trend calculation function."""
    trends = calculate_trend(sample_mutation_data, sample_metadata, p_threshold=0.05)
    
    # Check basic structure
    assert isinstance(trends, pd.DataFrame)
    assert 'mutation_id' in trends.columns
    assert 'slope' in trends.columns
    assert 'intercept' in trends.columns
    assert 'r_squared' in trends.columns
    assert 'p_value' in trends.columns
    assert 'significant' in trends.columns
    
    # Check trend calculations
    mut1_trend = trends[trends['mutation_id'] == 'mut1'].iloc[0]
    assert mut1_trend['slope'] > 0  # Increasing trend
    assert mut1_trend['significant'] == True
    
    mut2_trend = trends[trends['mutation_id'] == 'mut2'].iloc[0]
    assert mut2_trend['slope'] < 0  # Decreasing trend
    assert mut2_trend['significant'] == True

def test_calculate_trend_empty_input():
    """Test trend calculation with empty input."""
    empty_mutations = pd.DataFrame(columns=['mutation_id', 'sample_id', 'frequency'])
    empty_metadata = pd.DataFrame(columns=['sample_id', 'timepoint'])
    
    trends = calculate_trend(empty_mutations, empty_metadata, p_threshold=0.05)
    assert len(trends) == 0
    assert all(col in trends.columns for col in [
        'mutation_id', 'slope', 'intercept', 'r_squared',
        'p_value', 'significant'
    ])

def test_prepare_trajectory(sample_mutation_data, sample_metadata):
    """Test trajectory preparation function."""
    trajectories = prepare_trajectory(sample_mutation_data, sample_metadata)
    
    # Check basic structure
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
    
    # Check trajectory calculations
    mut1_traj = trajectories[trajectories['mutation_id'] == 'mut1'].iloc[0]
    assert mut1_traj['slope'] > 0
    assert mut1_traj['trajectory_type'] == 'increasing'
    assert mut1_traj['frequency_change'] > 0
    
    mut2_traj = trajectories[trajectories['mutation_id'] == 'mut2'].iloc[0]
    assert mut2_traj['slope'] < 0
    assert mut2_traj['trajectory_type'] == 'decreasing'
    assert mut2_traj['frequency_change'] < 0

def test_prepare_trajectory_empty_input():
    """Test trajectory preparation with empty input."""
    empty_mutations = pd.DataFrame(columns=['mutation_id', 'sample_id', 'frequency'])
    empty_metadata = pd.DataFrame(columns=['sample_id', 'timepoint'])
    
    trajectories = prepare_trajectory(empty_mutations, empty_metadata)
    assert len(trajectories) == 0
    assert all(col in trajectories.columns for col in [
        'mutation_id', 'slope', 'intercept', 'r_squared',
        'p_value', 'mean_frequency', 'std_frequency',
        'frequency_change', 'trajectory_type'
    ]) 