#!/usr/bin/env python3
"""
Filter mutations based on various criteria.

This module provides functions to filter mutations based on coverage.
It processes mutation data from inStrain output and applies filtering
criteria.

Inputs
------
* ``snv_info.tsv``      : Raw mutation data from inStrain
* ``processed_scaffolds.tsv`` : Filtered scaffolds that passed QC
* ``metadata.txt``      : Sample-level metadata

Outputs
-------
* ``filtered_mutations.tsv`` : Mutations that pass all filters
"""

# ─── Imports ──────────────────────────────────────────────────────────────────
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import pandas as pd

from utils import (
    setup_logging,
    get_logger,
    PerformanceMonitor,
)

# ─── Initialise a *temporary* logger ──────────────────────────────────────────
# It will start writing *after* proper configuration in `main()`.
logger = get_logger(__name__)
perf_monitor = PerformanceMonitor(logger)

# ─── Helper functions --------------------------------------------------------
def load_mutation_data(snv_file: Path | str) -> pd.DataFrame:
    """Load mutation data from a TSV file.

    Parameters
    ----------
    snv_file : Path
        Path to TSV file containing mutation data with columns

        * scaffold (str)            – scaffold identifier
        * position (int)            – 1-based coordinate
        * ref_base / new_base (str) – nucleotide bases
        * position_coverage (int)   – read coverage
        * frequency (float)         – mutation frequency
        * Sample (str)              – sample identifier

    Returns
    -------
    pd.DataFrame
        Mutation table.
    """
    snv_file = Path(snv_file)
    logger.info("Loading mutation data from %s", snv_file)
    df = pd.read_csv(snv_file, sep="\t")
    if "Sample" not in df.columns:
        df["Sample"] = snv_file.parent.parent.name
    return df


def filter_by_coverage(
    mutations: pd.DataFrame, min_coverage: int = 10
) -> pd.DataFrame:
    """Keep rows whose coverage is ≥ ``min_coverage``."""
    logger.info("Filtering mutations by minimum coverage of %s×", min_coverage)
    return mutations.loc[mutations["position_coverage"] >= min_coverage]


def load_metadata(meta_csv: Path, *, log: Optional[logging.Logger] = None) -> pd.DataFrame:
    """Load and pre-clean metadata from a CSV file.

    Raises
    ------
    ValueError
        If any required columns are missing.
    """
    # Required / optional fields ------------------------------------------------
    essential_cols = [
        "External.ID",
        "week_num",
        "Participant ID",
        "diagnosis",
    ]
    optional_cols = [
        "sex",
        "Height",
        "Weight",
        "BMI",
        "fecalcal_ng_ml",
        "Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)",
        "Antibiotics",
        "Immunosuppressants (e.g. oral corticosteroids)",
        "hbi",
        "sccai",
        "calprotectin_less_than_50", "blood_crp_less_than_1", "lysozyme_less_than_600", "secretory_iga_50_to_200"
    ]

    # Inspect header first so we can trim eagerly --------------------------------
    with open(meta_csv, "r") as f:
        first_line = f.readline().strip()
        # Detect delimiter: tab or comma
        if '\t' in first_line:
            delimiter = '\t'
            header = first_line.split('\t')
        else:
            delimiter = ','
            header = first_line.split(',')
    
    meta_keep = [c for c in essential_cols if c in header]
    missing_essentials = [c for c in essential_cols if c not in header]
    if missing_essentials:
        raise ValueError(f"Missing required metadata columns: {missing_essentials}")

    meta_keep += [c for c in optional_cols if c in header]
    missing_optional = [c for c in optional_cols if c not in header]
    (log or logger).warning(
        "Optional metadata columns missing and will be skipped: %s",
        missing_optional,
    )

    meta = pd.read_csv(meta_csv, usecols=meta_keep, dtype=str, sep=delimiter)
    # Normalise column names
    meta = meta.rename(
        columns={
            "Alcohol": "Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)",
            "Immunosuppressants": "Immunosuppressants (e.g. oral corticosteroids)",
        }
    )
    return meta


def filter_mutations(
    snv_file: Path,
    output_file: Path,
    metadata_file: Path,
    processed_scaffolds_file: Path,
    min_coverage: int = 10,
    min_week_occurrence: int = 3,
) -> None:
    """Core pipeline for filtering, merging, and writing mutations.
    
    The pipeline applies the following filters:
    1. Coverage filter: Keep mutations with coverage >= min_coverage
    2. Scaffold filter: Keep mutations only on processed scaffolds
    3. Week occurrence filter: Keep mutations that occur in >= min_week_occurrence weeks
    4. Duplicate removal: Remove duplicate scaffold+position+Sample combinations
    """
    # 1) Load -------------------------------------------------------------------
    mutations = load_mutation_data(snv_file)
    logger.info("Loaded %s rows from SNV file", len(mutations))

    # 2) Filter by coverage -----------------------------------------------------
    mutations = filter_by_coverage(mutations, min_coverage)
    logger.info("%s rows remain after coverage filter", len(mutations))

    # 3) Merge with processed scaffolds ----------------------------------------
    processed_scaffolds = pd.read_csv(processed_scaffolds_file, sep="\t")
    mutations = mutations.merge(
        processed_scaffolds, on=["scaffold", "Sample"], how="inner"
    )
    logger.info("%s rows remain after scaffold merge", len(mutations))

    # 4) Merge with metadata ----------------------------------------------------
    meta = load_metadata(metadata_file)
    logger.info(
        "Merging mutations (%s distinct Sample) with metadata (%s distinct External.ID)",
        mutations["Sample"].nunique(),
        meta["External.ID"].nunique(),
    )
    mutations = mutations.merge(
        meta, left_on="Sample", right_on="External.ID", how="inner"
    )
    logger.info("%s rows remain after metadata merge", len(mutations))

    # 5) Filter by minimum occurrence across weeks -----------------------------
    # This filter ensures mutations are only kept if they occur in at least 
    # min_week_occurrence different weeks for the same scaffold+position.
    # This helps filter out sequencing errors, rare mutations, and noise.
    # Count distinct weeks per scaffold+position combination
    week_counts = mutations.groupby(['scaffold', 'position'])['week_num'].nunique()
    mutations = mutations.merge(
        week_counts.reset_index().rename(columns={'week_num': 'week_count'}),
        on=['scaffold', 'position'],
        how='inner'
    )
    
    # Filter to keep only mutations that occur in at least min_week_occurrence weeks
    before = len(mutations)
    mutations = mutations[mutations['week_count'] >= min_week_occurrence]
    logger.info(
        "Filtered to mutations occurring in ≥%s weeks: %s rows remain (dropped %s)",
        min_week_occurrence, len(mutations), before - len(mutations)
    )
    
    # Drop the temporary week_count column
    mutations = mutations.drop(columns=['week_count'])
    
    # 6) Drop duplicates --------------------------------------------------------
    before = len(mutations)
    mutations = mutations.drop_duplicates(
        subset=["scaffold", "position", "Sample"], keep="first"
    )
    logger.info("Dropped %s duplicate rows", before - len(mutations))

    # 6) Write ------------------------------------------------------------------
    logger.info("Writing %s filtered mutations to %s", len(mutations), output_file)
    mutations.to_csv(output_file, sep="\t", index=False)


# ─── CLI / main ───────────────────────────────────────────────────────────────


def parse_cli():
    """Parse command-line arguments and return a Namespace."""
    import argparse

    p = argparse.ArgumentParser(
        description="Filter mutations based on coverage and sample metadata"
    )
    p.add_argument("--snv-file", type=Path, required=True, help="Input SNV .tsv")
    p.add_argument(
        "--output-file", type=Path, required=True, help="Filtered mutation output"
    )
    p.add_argument("--metadata-file", type=Path, required=True, help="Metadata .csv")
    p.add_argument(
        "--processed-scaffolds-file",
        type=Path,
        required=True,
        help="Processed scaffolds .tsv",
    )
    p.add_argument(
        "--min-coverage", type=int, default=10, help="Minimum coverage threshold"
    )
    p.add_argument(
        "--min-week-occurrence", type=int, default=3, 
        help="Minimum number of weeks a mutation must occur in (default: 3)"
    )
    p.add_argument(
        "--log-file",
        type=Path,
        help="Write detailed log to this file (and still echo to stderr)",
    )
    return p.parse_args()


def main() -> None:
    """Entry-point when executed as a script."""
    args = parse_cli()

    # ─── Proper logging config *before* grabbing module logger ────────────────
    setup_logging(args.log_file)
    global logger  # refresh the module-level logger so it inherits handlers
    logger = logging.getLogger(__name__)
    logger.info("Logging initialised - writing to %s", args.log_file or "stderr")

    # ─── Run filtering ────────────────────────────────────────────────────────
    perf_monitor.start_operation("filter_mutations")
    filter_mutations(
        snv_file=args.snv_file,
        output_file=args.output_file,
        metadata_file=args.metadata_file,
        processed_scaffolds_file=args.processed_scaffolds_file,
        min_coverage=args.min_coverage,
        min_week_occurrence=args.min_week_occurrence,
    )
    perf_monitor.end_operation("filter_mutations")


# ─── Execute ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
