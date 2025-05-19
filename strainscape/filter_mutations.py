#!/usr/bin/env python3
"""
Filter mutations based on various criteria.

This module provides functions to filter mutations based on coverage. 
It processes mutation data from
inStrain output and applies filtering criteria.

Inputs:
  - snv_info.tsv: Raw mutation data from inStrain
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
    PerformanceMonitor
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
    return mutations[mutations['position_coverage'] >= min_coverage]

def load_metadata(meta_csv, log=None):
    """Load and process metadata from CSV file.
    
    Args:
        meta_csv: Path to metadata CSV file
        log: Logger instance
        
    Returns:
        DataFrame containing processed metadata
    """
    # Define metadata columns to keep
    essential_cols = [
        'External.ID',
        'week_num',
        'Participant ID',
        'diagnosis'
    ]
    optional_cols = [
        'sex',
        'Height',
        'Weight',
        'BMI',
        'fecalcal_ng_ml',
        'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
        'Antibiotics',
        'Immunosuppressants (e.g. oral corticosteroids)',
        'hbi',
        'sccai'
    ]
    # Read the header to determine which columns are present
    with open(meta_csv, 'r') as f:
        header = f.readline().strip().split(',')
    meta_keep = [col for col in essential_cols if col in header]
    missing_essentials = [col for col in essential_cols if col not in header]
    if missing_essentials:
        raise ValueError(f"Missing required metadata columns: {missing_essentials}")
    # Add optional columns that are present
    meta_keep += [col for col in optional_cols if col in header]
    missing_optional = [col for col in optional_cols if col not in header]
    if missing_optional:
        log = log or logger
        log.warning(f"Optional metadata columns missing and will be skipped: {missing_optional}")
    # Load metadata
    meta = pd.read_csv(meta_csv, usecols=meta_keep, dtype=str)
    
    # Rename columns to match test data
    meta = meta.rename(columns={
        'Alcohol': 'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
        'Immunosuppressants': 'Immunosuppressants (e.g. oral corticosteroids)'
    })
    return meta

def filter_mutations(snv_file: Path,
                    output_file: Path,
                    metadata_file: Path,
                    processed_scaffolds_file: Path,
                    min_coverage: int = 10) -> None:
    """Filter mutations based on coverage and frequency criteria.
    
    Args:
        snv_file: Path to TSV file containing mutation data and coverage data
        output_file: Path to write filtered mutations
        metadata_file: Path to metadata CSV file
        processed_scaffolds_file: Path to processed scaffolds file
        min_coverage: Minimum required coverage (default: 10)
        
    Returns:
        None. Writes filtered mutations to output_file.
    """
    # Load data
    mutations = load_mutation_data(snv_file)
    logger.info(f"Reading {len(mutations):,}")
    
    # Apply filters
    mutations = filter_by_coverage(mutations, min_coverage)
    logger.info(f"Number of mutations after filtering: {len(mutations):,}")

    # Load processed scaffolds
    processed_scaffolds = pd.read_csv(processed_scaffolds_file, sep='\t')
    # Merge with processed scaffolds using inner join on both scaffold and Sample
    mutations = mutations.merge(processed_scaffolds, on=["scaffold", "Sample"], how="inner")
    logger.info(f"Number of mutations after merging: {len(mutations):,}")
    
    # Load and merge metadata
    meta = load_metadata(metadata_file)
    
    # Check for duplicates before merge
    logger.info(f"Number of unique Sample values in mutations: {mutations['Sample'].nunique():,}")
    logger.info(f"Number of unique External.ID values in metadata: {meta['External.ID'].nunique():,}")
    
    # Merge with metadata
    mutations = mutations.merge(meta, left_on="Sample", right_on="External.ID", how="inner")
    logger.info(f"Number of mutations after merging: {len(mutations):,}")
    
    # Drop duplicates keeping the first occurrence
    mutations = mutations.drop_duplicates(subset=['scaffold', 'position', 'Sample'], keep='first')
    logger.info(f"Number of mutations after dropping duplicates: {len(mutations):,}")
    
    # Write output
    logger.info(f"Writing {len(mutations):,} filtered mutations to {output_file}")
    mutations.to_csv(output_file, sep='\t', index=False)

def main():
    """Main function to filter mutations."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter mutations based on coverage and frequency')
    parser.add_argument('--snv-file', type=Path, required=True,
                      help='Path to SNV information TSV file')
    parser.add_argument('--output-file', type=Path, required=True,
                      help='Path to write filtered mutations')
    parser.add_argument('--metadata-file', type=Path, required=True,
                      help='Path to metadata CSV file')
    parser.add_argument('--processed-scaffolds-file', type=Path, required=True,
                      help='Path to processed scaffolds file')
    parser.add_argument('--min-coverage', type=int, default=10,
                      help='Minimum required coverage')
    parser.add_argument('--log-file', type=Path,
                      help='Path to write log file')
    
    args = parser.parse_args()
    
    # Setup logging
    if args.log_file:
        setup_logging(args.log_file)
    
    # Filter mutations
    filter_mutations(
        args.snv_file,
        args.output_file,
        args.metadata_file,
        args.processed_scaffolds_file,
        args.min_coverage
    )

if __name__ == '__main__':
    main()
