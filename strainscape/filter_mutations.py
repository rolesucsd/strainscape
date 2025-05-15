#!/usr/bin/env python3
"""
Filter mutations based on various criteria.

This module provides functions to filter mutations based on coverage,
frequency, and other quality metrics. It processes mutation data from
inStrain output and applies filtering criteria to identify high-quality
mutations for downstream analysis.

Inputs:
  - snv_info.tsv: Raw mutation data from inStrain
  - coverage.tsv: Coverage information
Outputs:
  - filtered_mutations.tsv: Filtered mutation data
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import logging

from strainscape.utils import (
    setup_logging,
    get_logger,
    PerformanceMonitor,
    validate_file_exists,
    read_large_csv
)

# Set up logging
logger = get_logger(__name__)
perf_monitor = PerformanceMonitor(logger)

def load_mutation_data(snv_file: Path) -> pd.DataFrame:
    """Load mutation data from a TSV file.
    
    Args:
        snv_file: Path to TSV file containing mutation data with columns:
            - scaffold: str, scaffold identifier
            - position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
            - coverage: int, read coverage
            - frequency: float, mutation frequency
            
    Returns:
        DataFrame containing mutation data
    """
    logger.info(f"Loading mutation data from {snv_file}")
    return pd.read_csv(snv_file, sep='\t')

def load_coverage_data(coverage_file: Path) -> pd.DataFrame:
    """Load coverage data from a TSV file.
    
    Args:
        coverage_file: Path to TSV file containing coverage data with columns:
            - scaffold: str, scaffold identifier
            - position: int, 1-based position
            - coverage: int, read coverage
            
    Returns:
        DataFrame containing coverage data
    """
    logger.info(f"Loading coverage data from {coverage_file}")
    return pd.read_csv(coverage_file, sep='\t')

def filter_by_coverage(mutations: pd.DataFrame,
                      min_coverage: int = 10) -> pd.DataFrame:
    """Filter mutations by minimum coverage threshold.
    
    Args:
        mutations: DataFrame containing mutation data
        min_coverage: Minimum required coverage (default: 10)
        
    Returns:
        DataFrame containing mutations that pass coverage filter
    """
    logger.info(f"Filtering mutations by minimum coverage of {min_coverage}")
    return mutations[mutations['coverage'] >= min_coverage]

def filter_by_frequency(mutations: pd.DataFrame,
                       min_freq: float = 0.05,
                       max_freq: float = 0.95) -> pd.DataFrame:
    """Filter mutations by frequency range.
    
    Args:
        mutations: DataFrame containing mutation data
        min_freq: Minimum allowed frequency (default: 0.05)
        max_freq: Maximum allowed frequency (default: 0.95)
        
    Returns:
        DataFrame containing mutations that pass frequency filter
    """
    logger.info(f"Filtering mutations by frequency range [{min_freq}, {max_freq}]")
    return mutations[
        (mutations['frequency'] >= min_freq) &
        (mutations['frequency'] <= max_freq)
    ]

def filter_mutations(snv_file: Path,
                    coverage_file: Path,
                    output_file: Path,
                    min_coverage: int = 10,
                    min_freq: float = 0.05,
                    max_freq: float = 0.95) -> None:
    """Filter mutations based on coverage and frequency criteria.
    
    Args:
        snv_file: Path to TSV file containing mutation data
        coverage_file: Path to TSV file containing coverage data
        output_file: Path to write filtered mutations
        min_coverage: Minimum required coverage (default: 10)
        min_freq: Minimum allowed frequency (default: 0.05)
        max_freq: Maximum allowed frequency (default: 0.95)
        
    Returns:
        None. Writes filtered mutations to output_file.
    """
    # Load data
    mutations = load_mutation_data(snv_file)
    coverage = load_coverage_data(coverage_file)
    
    # Merge with coverage data
    mutations = mutations.merge(coverage, on=['scaffold', 'position'], how='left')
    
    # Apply filters
    mutations = filter_by_coverage(mutations, min_coverage)
    mutations = filter_by_frequency(mutations, min_freq, max_freq)
    
    # Write output
    logger.info(f"Writing {len(mutations):,} filtered mutations to {output_file}")
    mutations.to_csv(output_file, sep='\t', index=False)

def filter_mutations_df(mutations: pd.DataFrame, min_coverage: float = 10, min_freq: float = 0.0) -> pd.DataFrame:
    """Filter mutations DataFrame by coverage and frequency (for unit tests)."""
    return mutations[(mutations['coverage'] >= min_coverage) & (mutations['frequency'] >= min_freq)].copy()

def main():
    """Main function to filter mutations."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter mutations based on coverage and frequency')
    parser.add_argument('--snv-file', type=Path, required=True,
                      help='Path to SNV information TSV file')
    parser.add_argument('--coverage-file', type=Path, required=True,
                      help='Path to coverage information TSV file')
    parser.add_argument('--output-file', type=Path, required=True,
                      help='Path to write filtered mutations')
    parser.add_argument('--min-coverage', type=int, default=10,
                      help='Minimum required coverage')
    parser.add_argument('--min-freq', type=float, default=0.05,
                      help='Minimum allowed frequency')
    parser.add_argument('--max-freq', type=float, default=0.95,
                      help='Maximum allowed frequency')
    parser.add_argument('--log-file', type=Path,
                      help='Path to write log file')
    
    args = parser.parse_args()
    
    # Setup logging
    if args.log_file:
        setup_logging(args.log_file)
    
    # Filter mutations
    filter_mutations(
        args.snv_file,
        args.coverage_file,
        args.output_file,
        args.min_coverage,
        args.min_freq,
        args.max_freq
    )

if __name__ == '__main__':
    main() 