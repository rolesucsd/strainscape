#!/usr/bin/env python3
"""
Process scaffold information and map mutations to scaffolds.

This module provides functions to process scaffold information and map mutations
to their corresponding scaffolds. It handles scaffold mapping files and
generates detailed output about mutation locations relative to scaffolds.

Inputs:
  - scaffold_info.tsv: Scaffold mapping information
  - snv_info.tsv: Mutation information
  - assembly.stb: Scaffold to bin mapping
Outputs:
  - mapped_mutations.tsv: Mutations mapped to scaffolds with additional information
"""

from pathlib import Path
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse, logging, sys, textwrap
from typing import Dict, List, Optional, Tuple, Union
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

# ────────── main pipeline ──────────
def process(scaffold_info, bin_file, output_file, min_length=1000,
            min_coverage=5.0, min_breadth=0.4,
            min_completeness=50.0, max_contamination=10.0,
            log=None):
    """Process scaffolds with quality filters.
    
    Args:
        scaffold_info: Path to scaffold info file
        bin_file: Path to bin.txt file (scaffold-bin mapping)
        output_file: Path to output file
        min_length: Minimum scaffold length
        min_coverage: Minimum coverage
        min_breadth: Minimum breadth
        log: Logger instance
    """
    if log is None:                         # Fallback logger
        log = logger

    # 1. bin map -------------------------------------------------------------
    bin_df = pd.read_csv(
        bin_file,
        sep="\t",
        dtype={"scaffold": str, "bin": str, "Completeness": float, "Contamination": float, "Genome_Size": float}
    )

    log.info(f"Loaded {len(bin_df):,} scaffold↦bin mappings from {bin_file}")

    # 3. scaffold_info ------------------------------------------------------
    desired_cols = ["scaffold", "length", "coverage", "breadth", "nucl_diversity", "Sample"]
    with open(scaffold_info, "r") as f:
        header = f.readline().strip().split("\t")
    usecols = [c for c in desired_cols if c in header]
    dtypes = {
        "scaffold": str,
        "length": np.int32,
        "coverage": np.float32,
        "breadth": np.float32,
    }
    scaff = pd.read_csv(scaffold_info, sep="\t", usecols=usecols, dtype=dtypes)
    scaff = scaff.merge(bin_df, on="scaffold", how="left")

    log.info(f"Scaffolds loaded: {len(scaff):,}")
    query = "length >= @min_length and coverage >= @min_coverage and breadth >= @min_breadth"
    if "Completeness" in scaff.columns and "Contamination" in scaff.columns:
        query += " and Completeness >= @min_completeness and Contamination <= @max_contamination"
    scaff = scaff.query(query)
    log.info(f"After quality filters: {len(scaff):,}")

# 6. summary ------------------------------------------------------------
    summary_cols = ["length", "coverage", "breadth"]
    for col in ["Completeness", "Contamination", "Genome_Size"]:
        if col in scaff.columns:
            summary_cols.append(col)
    summary = (
        scaff.groupby("bin")[summary_cols]
        .agg(["mean", "std", "min", "max"])
        .round(2)
    )
    log.info("Summary per bin:\n%s", summary)

    # 7. write --------------------------------------------------------------
    scaff.to_csv(output_file, sep="\t", index=False)
    log.info(f"Written {len(scaff):,} rows → {output_file}")


def process_scaffolds(
    scaffold_file: str,
    bin_file: str | None = None,
    output_file: str | None = None,
    min_length: int = 1000,
    min_coverage: float = 5.0,
    min_breadth: float = 0.4,
    min_completeness: float = 50,
    max_contamination: float = 10,
    log_file: Optional[str] = None,
    *,
    metadata_file: str | None = None,
    bin_dir: str | None = None,
    threads: int = 1,
    chunksize: int = 10000,
    **_: object,
) -> None:
    """
    Process scaffolds with quality filters.
    
    Args:
        scaffold_file: Path to scaffold info TSV
        bin_file: Path to bin.txt file
        output_file: Path to output file
        min_length: Minimum scaffold length
        min_coverage: Minimum coverage
        min_breadth: Minimum breadth
        log_file: Optional path to log file
    """
    if log_file:
        setup_logging(log_file)

    if bin_file is None:
        raise ValueError("bin_file or bin_dir must be provided")

    scaffold_path = Path(scaffold_file)

    process(
        scaffold_info=str(scaffold_path),
        bin_file=bin_file,
        output_file=output_file,
        min_length=min_length,
        min_coverage=min_coverage,
        min_breadth=min_breadth,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        log=logger
    )

    if metadata_file:
        df = pd.read_csv(output_file, sep='\t')
        meta = pd.read_csv(metadata_file, dtype=str)
        df = df.merge(meta, left_on='Sample', right_on='External.ID', how='left')
        df.to_csv(output_file, sep='\t', index=False)

    if bin_dir:
        os.remove(bin_file)


# ────────── CLI ──────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent("""
            Merge inStrain scaffold_info with STB, bin.txt
            Provides the same output as the original script but with less I/O
            and no intermediate files.
        """))
    ap.add_argument("--scaffold_file", required=True)
    ap.add_argument("--bin_file",      required=True)
    ap.add_argument("--output_file",   required=True)
    ap.add_argument("--min_length",   type=int,   default=1000)
    ap.add_argument("--min_coverage", type=float, default=5.0)
    ap.add_argument("--min_breadth",  type=float, default=0.4)
    ap.add_argument("--min_completeness", type=float, default=50.0)
    ap.add_argument("--max_contamination", type=float, default=10.0)
    ap.add_argument("--log_file",     required=False)
    args = ap.parse_args()

    process_scaffolds(
        scaffold_file=args.scaffold_file,
        bin_file=args.bin_file,
        output_file=args.output_file,
        min_length=args.min_length,
        min_coverage=args.min_coverage,
        min_breadth=args.min_breadth,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        log_file=args.log_file
    )
