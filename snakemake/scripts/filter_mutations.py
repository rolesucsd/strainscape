#!/usr/bin/env python3
"""
Merge SNV information with STB, metadata, and filtered scaffolds.

This script combines SNV information with bin mapping (STB), patient metadata,
and filtered scaffold information to ensure SNVs only map to scaffolds that
passed previous filtering criteria.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import logging

from utils import (
    setup_logging,
    PerformanceMonitor,
    validate_file_exists,
    read_large_csv
)

# Set up logging
logger = setup_logging()
perf_monitor = PerformanceMonitor(logger)

def merge_snv_info(
    snv_file: str,
    processed_scaffold_file: str,
    output_file: str,
    threads: int = 4,
    chunksize: int = 10000
) -> None:
    """
    Merge SNV information with processed scaffolds on 'scaffold'.
    Only keep SNVs mapping to scaffolds present in processed_scaffolds.tsv.
    """
    perf_monitor.start_operation('merge_snv_info')
    logger.info(f"Processing files: {snv_file}, {processed_scaffold_file}")
    try:
        validate_file_exists(snv_file)
        validate_file_exists(processed_scaffold_file)

        # Load processed scaffold information
        logger.info("Loading processed scaffold information")
        processed_scaffolds = pd.read_csv(processed_scaffold_file, sep='\t')
        valid_scaffolds = set(processed_scaffolds['scaffold'].unique())
        logger.info(f"Found {len(valid_scaffolds)} valid scaffolds")

        # Process SNV information in chunks
        chunks = []
        for chunk in pd.read_csv(snv_file, chunksize=chunksize, sep='\t'):
            chunk.columns = chunk.columns.str.strip()
            # Keep only SNVs mapping to valid scaffolds and with coverage >= 10
            chunk = chunk[(chunk['scaffold'].isin(valid_scaffolds)) & (chunk['position_coverage'] >= 10)]
            if not chunk.empty:
                chunks.append(chunk)

        # Combine and save results
        if chunks:
            merged_snvs = pd.concat(chunks, ignore_index=True)
            merged_snvs.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved {len(merged_snvs)} filtered SNVs to {output_file}")
            logger.info(f"Number of unique scaffolds: {merged_snvs['scaffold'].nunique()}")
        else:
            logger.warning("No SNV information found after filtering")
            pd.DataFrame().to_csv(output_file, sep='\t', index=False)
        perf_monitor.end_operation('merge_snv_info')
    except Exception as e:
        logger.error(f"Error merging SNV information: {str(e)}")
        raise

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Filter SNVs to only those mapping to scaffolds present in processed_scaffolds.tsv.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--snv_file", required=True, help="Path to SNV information file")
    parser.add_argument("--processed_scaffold_file", required=True, help="Path to processed scaffold file")
    parser.add_argument("--output_file", required=True, help="Path to output file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--chunksize", type=int, default=10000, help="Size of chunks for processing")
    parser.add_argument("--log_file", help="Path to log file (optional)")
    args = parser.parse_args()
    if args.log_file:
        logger = setup_logging(args.log_file)
    try:
        merge_snv_info(
            args.snv_file,
            args.processed_scaffold_file,
            args.output_file,
            args.threads,
            args.chunksize
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise 