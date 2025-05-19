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
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse, logging, sys, textwrap
from typing import Dict, List, Optional, Tuple, Union
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

# ────────── main pipeline ──────────
def process(scaffold_info, bin_file, output_file, min_length=1000, min_coverage=5.0, min_breadth=0.4, log=None):
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
    bin_df = pd.read_csv(bin_file, sep="\t", names=["scaffold", "bin"], dtype={"scaffold": str, "bin": str})
    log.info(f"Loaded {len(bin_df):,} scaffold↦bin mappings from {bin_file}")

    # 3. scaffold_info ------------------------------------------------------
    usecols = ["scaffold", "length", "coverage", "breadth", "nucl_diversity", "Sample"]
    dtypes  = {"scaffold": str,
               "length":   np.int32,
               "coverage": np.float32,
               "breadth":  np.float32}
    scaff = pd.read_csv(scaffold_info, sep="\t", usecols=usecols, dtype=dtypes)
    scaff = scaff.merge(bin_df, on="scaffold", how="left")

    log.info(f"Scaffolds loaded: {len(scaff):,}")
    scaff = scaff.query("length >= @min_length and coverage >= @min_coverage and breadth >= @min_breadth")
    log.info(f"After quality filters: {len(scaff):,}")

    # 6. summary ------------------------------------------------------------
    summary = (scaff
               .groupby("bin")[["length", "coverage", "breadth"]]
               .agg(["mean", "std", "min", "max"])
               .round(2))
    log.info("Summary per bin:\n%s", summary)

    # 7. write --------------------------------------------------------------
    scaff.to_csv(output_file, sep="\t", index=False)
    log.info(f"Written {len(scaff):,} rows → {output_file}")


def process_scaffolds(
    scaffold_file: str,
    bin_file: str,
    output_file: str,
    min_length: int = 1000,
    min_coverage: float = 5.0,
    min_breadth: float = 0.4,
    log_file: Optional[str] = None
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
    
    process(
        scaffold_info=scaffold_file,
        bin_file=bin_file,
        output_file=output_file,
        min_length=min_length,
        min_coverage=min_coverage,
        min_breadth=min_breadth,
        log=logger
    )


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
    ap.add_argument("--log_file",     required=True)
    args = ap.parse_args()

    process_scaffolds(
        scaffold_file=args.scaffold_file,
        bin_file=args.bin_file,
        output_file=args.output_file,
        min_length=args.min_length,
        min_coverage=args.min_coverage,
        min_breadth=args.min_breadth,
        log_file=args.log_file
    )
