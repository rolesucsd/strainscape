#!/usr/bin/env python3
"""
Generate a bin.txt file from a bin folder.

This script reads all FASTA files in a bin folder and writes a tab-separated
bin.txt file with columns: scaffold, bin.

Usage:
    python make_bin_txt.py --bin_dir /path/to/bin/folder --output_file /path/to/output/bin.txt
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import argparse, logging, sys, textwrap
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

def make_bin_txt(bin_dir: Path, output_file: str, log: logging.Logger = None) -> None:
    """
    Generate a bin.txt file from a bin folder.

    Args:
        bin_dir: Path to the bin folder containing FASTA files
        output_file: Path to the output bin.txt file
        log: Logger instance
    """
    if log is None:
        log = logger

    rows = []
    for fa in bin_dir.glob("*.fa"):
        bin_name = fa.stem
        for record in SeqIO.parse(fa, "fasta"):
            rows.append((record.id, bin_name))
    df = pd.DataFrame(rows, columns=["scaffold", "bin"])
    log.info(f"Collected {len(df):,} scaffold↦bin mappings from {bin_dir}")
    df.to_csv(output_file, sep="\t", index=False, header=False)
    log.info(f"Written {len(df):,} rows → {output_file}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent("""
            Generate a bin.txt file from a bin folder.
        """))
    ap.add_argument("--bin_dir", required=True, help="Path to the bin folder containing FASTA files")
    ap.add_argument("--output_file", required=True, help="Path to the output bin.txt file")
    ap.add_argument("--log_file", required=True, help="Path to the log file")
    args = ap.parse_args()

    setup_logging(args.log_file)
    make_bin_txt(Path(args.bin_dir), args.output_file, logger) 