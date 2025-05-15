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

# ────────── helper: build bin<->scaffold DF (no temp file) ──────────
def bin_dataframe(bin_dir: Path, log: logging.Logger) -> pd.DataFrame:
    rows = []
    for fa in bin_dir.glob("*.fa"):
        bin_name = fa.stem
        for record in SeqIO.parse(fa, "fasta"):
            rows.append((record.id, bin_name))
    df = pd.DataFrame(rows, columns=["scaffold", "bin"])
    log.info(f"Collected {len(df):,} scaffold↦bin mappings from {bin_dir}")
    return df


# ────────── main pipeline ──────────
def process(scaffold_info, stb_file, meta_csv, bin_dir, output_file, min_length=1000, min_coverage=5.0, min_breadth=0.4, log=None):
    """Process scaffolds with quality filters and metadata.
    
    Args:
        scaffold_info: Path to scaffold info file
        stb_file: Path to STB file
        meta_csv: Path to metadata CSV
        bin_dir: Directory containing bin files
        output_file: Path to output file
        min_length: Minimum scaffold length
        min_coverage: Minimum coverage
        min_breadth: Minimum breadth
        log: Logger instance
    """
    # Define metadata columns to keep
    meta_keep = [
        'External.ID',
        'week_num',
        'Participant ID',
        'sex',
        'diagnosis',
        'Height',
        'Weight',
        'BMI',
        'fecalcal_ng_ml',
        'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
        'Antibiotics',
        'Immunosuppressants (e.g. oral corticosteroids)'
    ]
    
    # Load metadata
    meta = pd.read_csv(meta_csv, usecols=meta_keep, dtype=str)
    
    # Rename columns to match test data
    meta = meta.rename(columns={
        'Alcohol': 'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
        'Immunosuppressants': 'Immunosuppressants (e.g. oral corticosteroids)'
    })

    if log is None:                         # Fallback logger
        log = logger

    # 1. bin map -------------------------------------------------------------
    bin_df = bin_dataframe(bin_dir, log)

    # 2. STB ----------------------------------------------------------------
    stb_df = pd.read_csv(stb_file, sep="\t",
                         names=["scaffold", "Sample"],  # Sample == External.ID
                         dtype={"scaffold": str, "Sample": str})

    stb_df = stb_df.merge(bin_df, on="scaffold", how="left")

    # 3. scaffold_info ------------------------------------------------------
    usecols = ["scaffold", "length", "coverage", "breadth"]
    dtypes  = {"scaffold": str,
               "length":   np.int32,
               "coverage": np.float32,
               "breadth":  np.float32}
    scaff = pd.read_csv(scaffold_info, sep="\t", usecols=usecols, dtype=dtypes)

    log.info(f"Scaffolds loaded: {len(scaff):,}")
    scaff = scaff.query("length >= @min_length and coverage >= @min_coverage and breadth >= @min_breadth")
    log.info(f"After quality filters: {len(scaff):,}")

    # 4. merge scaff <- stb+bin --------------------------------------------
    scaff = scaff.merge(stb_df, on="scaffold", how="inner")

    # 5. metadata  ----------------------------------------------------------
    scaff = scaff.merge(meta, left_on="Sample", right_on="External.ID", how="left")

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
    stb_file: str,
    metadata_file: str,
    bin_dir: str,
    output_file: str,
    min_length: int = 1000,
    min_coverage: float = 5.0,
    min_breadth: float = 0.4,
    threads: int = 4,
    chunksize: int = 10000,
    log_file: Optional[str] = None
) -> None:
    """
    Process scaffolds with quality filters and metadata.
    
    Args:
        scaffold_file: Path to scaffold info TSV
        stb_file: Path to STB file
        metadata_file: Path to metadata CSV
        bin_dir: Path to bin directory
        output_file: Path to output file
        min_length: Minimum scaffold length
        min_coverage: Minimum coverage
        min_breadth: Minimum breadth
        threads: Number of threads
        chunksize: Size of chunks for processing
        log_file: Optional path to log file
    """
    if log_file:
        setup_logging(log_file)
    
    process(
        scaffold_info=scaffold_file,
        stb_file=stb_file,
        meta_csv=metadata_file,
        bin_dir=bin_dir,
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
        description=textwrap.dedent("""\
            Merge inStrain scaffold_info with STB, bin FASTAs and sample metadata.
            Provides the same output as the original script but with less I/O
            and no intermediate files.
        """))
    ap.add_argument("--scaffold_file", required=True)
    ap.add_argument("--stb_file",      required=True)
    ap.add_argument("--metadata_file", required=True)
    ap.add_argument("--bin_dir",       required=True)
    ap.add_argument("--output_file",   required=True)
    ap.add_argument("--min_length",   type=int,   default=1000)
    ap.add_argument("--min_coverage", type=float, default=5.0)
    ap.add_argument("--min_breadth",  type=float, default=0.4)
    ap.add_argument("--log_file",     required=True)
    args = ap.parse_args()

    process_scaffolds(
        scaffold_file=args.scaffold_file,
        stb_file=args.stb_file,
        metadata_file=args.metadata_file,
        bin_dir=args.bin_dir,
        output_file=args.output_file,
        min_length=args.min_length,
        min_coverage=args.min_coverage,
        min_breadth=args.min_breadth,
        log_file=args.log_file
    )
