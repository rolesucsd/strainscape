#!/usr/bin/env python3
"""
Preprocess SNV file to be compatible with calculate_trends.py

This script converts SNV files with the format:
- sample, contig, position, ref_base, depth, A, C, G, T, ...

To the format expected by calculate_trends.py:
- scaffold, position, Sample, ref_base, position_coverage, week_num, A, C, G, T

It joins with metadata to add week_num information.
"""

from pathlib import Path
import argparse
import pandas as pd
from strainscape.utils import setup_logging, get_logger

logger = get_logger(__name__)


def preprocess_snvs(
    snv_file: Path,
    metadata_file: Path,
    output_file: Path,
    sample_col: str = "Sample",
    week_col: str = "week_num",
) -> None:
    """
    Preprocess SNV file to match calculate_trends.py format.

    Parameters
    ----------
    snv_file : Path
        Input SNV file with columns: sample, contig, position, ref_base, depth, A, C, G, T
    metadata_file : Path
        Metadata file with Sample and week_num columns
    output_file : Path
        Output file path for processed SNVs
    sample_col : str
        Column name in metadata that matches the sample column (default: "Sample")
    week_col : str
        Column name in metadata for week number (default: "week_num")
    """
    logger.info(f"Loading SNV file: {snv_file}")
    snvs = pd.read_csv(snv_file, sep="\t")

    # Check required columns exist
    required_cols = ["sample", "contig", "position", "ref_base", "depth", "A", "C", "G", "T"]
    missing_cols = [col for col in required_cols if col not in snvs.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns in SNV file: {missing_cols}. "
            f"Found columns: {list(snvs.columns)}"
        )

    logger.info(f"Loaded {len(snvs):,} rows from SNV file")

    # Rename columns to match expected format
    logger.info("Renaming columns...")
    snvs = snvs.rename(columns={
        "sample": "Sample",
        "contig": "scaffold",
        "depth": "position_coverage"
    })

    # Load metadata
    logger.info(f"Loading metadata file: {metadata_file}")
    metadata = pd.read_csv(metadata_file, sep=None, engine="python")
    
    # Check metadata has required columns
    if sample_col not in metadata.columns:
        raise ValueError(
            f"Metadata file missing '{sample_col}' column. "
            f"Found columns: {list(metadata.columns)}"
        )
    if week_col not in metadata.columns:
        raise ValueError(
            f"Metadata file missing '{week_col}' column. "
            f"Found columns: {list(metadata.columns)}"
        )

    logger.info(f"Loaded metadata with {len(metadata):,} rows")
    logger.info(f"Metadata columns: {list(metadata.columns)}")

    # Prepare metadata for merge (keep only necessary columns)
    meta_subset = metadata[[sample_col, week_col]].copy()
    
    # Convert week_num to numeric if it's not already
    meta_subset[week_col] = pd.to_numeric(meta_subset[week_col], errors="coerce")
    
    # Drop rows with missing week_num
    before_meta = len(meta_subset)
    meta_subset = meta_subset.dropna(subset=[week_col])
    if len(meta_subset) < before_meta:
        logger.warning(f"Dropped {before_meta - len(meta_subset)} metadata rows with missing {week_col}")

    # Merge with metadata
    logger.info(
        f"Merging SNVs ({snvs['Sample'].nunique()} distinct samples) "
        f"with metadata ({meta_subset[sample_col].nunique()} distinct samples)"
    )
    
    merged = snvs.merge(
        meta_subset,
        left_on="Sample",
        right_on=sample_col,
        how="inner"
    )
    
    logger.info(f"{len(merged):,} rows remain after metadata merge")

    # Drop rows without week_num (shouldn't happen with inner join, but just in case)
    before_week = len(merged)
    merged = merged.dropna(subset=[week_col])
    if len(merged) < before_week:
        logger.warning(f"Dropped {before_week - len(merged)} rows with missing {week_col}")

    # Rename week_col to week_num if needed
    if week_col != "week_num":
        merged = merged.rename(columns={week_col: "week_num"})

    # Drop the duplicate sample column from metadata if it exists
    if sample_col != "Sample" and sample_col in merged.columns:
        merged = merged.drop(columns=[sample_col])

    # Select only the columns needed by calculate_trends.py
    required_output_cols = [
        "scaffold", "position", "Sample", "ref_base", 
        "position_coverage", "week_num", "A", "C", "G", "T"
    ]
    
    # Keep any additional columns that were in the original file
    additional_cols = [col for col in merged.columns 
                      if col not in required_output_cols]
    
    output_cols = required_output_cols + additional_cols
    merged = merged[output_cols]

    # Write output
    logger.info(f"Writing {len(merged):,} rows to {output_file}")
    merged.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Successfully created preprocessed file: {output_file}")
    
    # Print summary statistics
    logger.info(f"Summary:")
    logger.info(f"  - Unique samples: {merged['Sample'].nunique()}")
    logger.info(f"  - Unique scaffolds: {merged['scaffold'].nunique()}")
    logger.info(f"  - Unique positions: {merged.groupby(['scaffold', 'position']).ngroups:,}")
    logger.info(f"  - Week range: {merged['week_num'].min():.0f} - {merged['week_num'].max():.0f}")


# ───────────── CLI ─────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Preprocess SNV file for calculate_trends.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("--snv_file", required=True,
                    help="Input SNV file (TSV with sample, contig, position, ref_base, depth, A, C, G, T)")
    ap.add_argument("--metadata_file", required=True,
                    help="Metadata file (CSV/TSV with Sample and week_num columns)")
    ap.add_argument("--output_file", required=True,
                    help="Output file path for preprocessed SNVs")
    ap.add_argument("--sample_col", default="Sample",
                    help="Column name in metadata that matches sample names")
    ap.add_argument("--week_col", default="week_num",
                    help="Column name in metadata for week number")
    ap.add_argument("--log_file", required=False,
                    help="Optional log file path")
    
    args = ap.parse_args()

    if args.log_file:
        setup_logging(args.log_file)

    preprocess_snvs(
        Path(args.snv_file),
        Path(args.metadata_file),
        Path(args.output_file),
        sample_col=args.sample_col,
        week_col=args.week_col,
    )

