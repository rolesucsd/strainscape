#!/usr/bin/env python3
"""
Generate a bin.txt file from a bin folder, optionally merging in CheckM2 metrics.

This script reads all FASTA files in a bin folder and writes a tab-separated
bin.txt file with columns: scaffold, bin, Completeness, Contamination, Genome_Size.

Usage:
    python make_bin_txt.py \
      --bin_dir /path/to/bin/folder \
      --checkm2-tsv /path/to/checkm2/quality_report.tsv \
      --output_file /path/to/output/bin.txt \
      --log_file /path/to/log.txt
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import argparse, logging, sys, textwrap
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

def make_bin_txt(
    bin_dir: Path,
    output_file: str,
    checkm2_tsv: Path | None = None,
    log: logging.Logger = None
) -> None:
    """
    Generate a bin.txt file from a bin folder and merge CheckM2 metrics if provided.

    Args:
        bin_dir: Path to the bin folder containing FASTA files
        output_file: Path to the output bin.txt file
        checkm2_tsv: Optional Path to CheckM2 quality_report.tsv
        log: Logger instance
    """
    if log is None:
        log = logger

    # 1) build scaffold-bin mapping
    rows = []
    for fa in bin_dir.glob("*.fa"):
        bin_name = fa.stem
        for record in SeqIO.parse(fa, "fasta"):
            rows.append((record.id, bin_name))
    df = pd.DataFrame(rows, columns=["scaffold", "bin"] )
    log.info(f"Collected {len(df):,} scaffold→bin mappings from {bin_dir}")

    # 2) merge with CheckM2 metrics if available
    if checkm2_tsv and checkm2_tsv.exists():
        log.info(f"Loading CheckM2 metrics from {checkm2_tsv}")
        quality_df = pd.read_csv(
            checkm2_tsv, sep="\t",
            usecols=["Name", "Completeness", "Contamination", "Genome_Size"]
        )
        merged = df.merge(
            quality_df,
            left_on="bin", right_on="Name",
            how="left"
        )
        merged = merged.drop(columns=["Name"] )
    else:
        merged = df

    # 3) write out merged file
    log.info(f"Writing {len(merged):,} rows → {output_file}")
    merged.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent(
            "Generate a bin.txt file from a bin folder, optionally merging CheckM2 metrics."
        )
    )
    ap.add_argument(
        "--bin_dir", required=True, type=Path,
        help="Path to the bin folder containing FASTA files"
    )
    ap.add_argument(
        "--checkm2-tsv", type=Path,
        help="Path to CheckM2 quality_report.tsv"
    )
    ap.add_argument(
        "--output_file", required=True,
        help="Path to the output bin.txt file"
    )
    ap.add_argument(
        "--log_file", required=True,
        help="Path to the log file"
    )
    args = ap.parse_args()

    setup_logging(args.log_file)
    make_bin_txt(
        args.bin_dir,
        args.output_file,
        args.checkm2_tsv,
        logger
    )
