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
import argparse, logging, sys, textwrap, gzip
from strainscape.utils import setup_logging, get_logger

logger = get_logger(__name__)

FA_EXTS   = ('.fa', '.fna', '.fasta')
GZ_EXTS   = tuple(e + '.gz' for e in FA_EXTS)
ALL_EXTS  = FA_EXTS + GZ_EXTS

def strip_ext(name: str) -> str:
    # remove the longest matching fasta/fna/fasta (+ optional .gz) suffix
    for ext in sorted(ALL_EXTS, key=len, reverse=True):
        if name.endswith(ext):
            return name[:-len(ext)]
    return Path(name).stem

def iter_fasta_files(bin_dir: Path):
    for p in bin_dir.iterdir():
        if any(p.name.endswith(ext) for ext in ALL_EXTS):
            yield p

def make_bin_txt(bin_dir: Path, output_file: str, checkm2_tsv: Path | None = None, log: logging.Logger = None) -> None:
    if log is None:
        log = logger

    rows = []
    files = list(iter_fasta_files(bin_dir))
    if not files:
        log.warning(f"No FASTA files found in {bin_dir} with extensions: {ALL_EXTS}")

    for fa in files:
        bin_name = strip_ext(fa.name)
        # robust opener for gz or plain
        opener = gzip.open if fa.suffix == '.gz' or fa.name.endswith('.gz') else open
        with opener(fa, 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                rows.append((record.id, bin_name))
    df = pd.DataFrame(rows, columns=["scaffold", "bin"])
    log.info(f"Collected {len(df):,} scaffold→bin mappings from {bin_dir}")

    if checkm2_tsv and checkm2_tsv.exists():
        log.info(f"Loading CheckM2 metrics from {checkm2_tsv}")
        quality_df = pd.read_csv(
            checkm2_tsv, sep="\t",
            usecols=["Name", "Completeness", "Contamination", "Genome_Size"]
        )
        # normalize CheckM2 Name the same way
        quality_df["Name"] = quality_df["Name"].apply(lambda s: strip_ext(Path(s).name))
        merged = df.merge(quality_df, left_on="bin", right_on="Name", how="left").drop(columns=["Name"])
    else:
        merged = df

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
