#!/usr/bin/env python3
"""
Calculate SNV bin statistics from all_snvs parquet files.

This script:
1. Reads all_snvs parquet files
2. Calculates SNVs per bin (patient_id, bin) combinations
3. Creates delta frequency histogram
4. Outputs summary statistics

Usage:
    python calculate_snv_bin_stats.py --input-dir <path> --output-dir <path> --meta-file <path>
"""

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(asctime)s %(message)s',
        datefmt='%H:%M:%S'
    )
    return logging.getLogger(__name__)

def load_metadata(meta_path: Path) -> pd.DataFrame:
    """Load and clean metadata."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading metadata from {meta_path}")
    
    meta = pd.read_csv(meta_path, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(' ', '_')
    
    # Handle different column names for patient ID
    if 'participant_id' in meta.columns and 'patient_id' not in meta.columns:
        meta = meta.rename(columns={'participant_id': 'patient_id'})
    
    # Keep only essential columns
    meta = meta[['patient_id', 'diagnosis']].drop_duplicates()
    meta['group'] = meta['diagnosis'].map({'nonIBD': 'nonIBD', 'UC': 'UC', 'CD': 'CD'})
    
    logger.info(f"Metadata loaded: {len(meta)} patients")
    logger.info(f"Group counts: {meta['group'].value_counts().to_dict()}")
    
    return meta

def load_snvs(input_dir: Path) -> pd.DataFrame:
    """Load all SNVs from parquet files."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading SNVs from {input_dir}")
    
    # Load all parquet files in the directory
    dataset = ds.dataset(input_dir, format='parquet')
    snvs = dataset.to_table().to_pandas()
    
    logger.info(f"SNVs loaded: {len(snvs):,} rows")
    logger.info(f"Columns: {list(snvs.columns)}")
    
    return snvs

def calculate_snvs_per_bin(snvs: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Calculate SNV counts per (patient_id, bin, group) in long form.

    Output columns: patient_id, bin, group, sweeps, non_sweeps, total
    """
    logger = logging.getLogger(__name__)
    logger.info("Calculating SNVs per bin...")

    # Merge SNVs with metadata to get group information
    snvs_with_meta = snvs.merge(meta, on='patient_id', how='left')

    # Count SNVs per bin, split by sweep status
    snv_counts = (
        snvs_with_meta
        .groupby(['patient_id', 'bin', 'group'])
        .agg({'is_sweep': ['sum', 'count']})  # sum = sweeps, count = total
        .round(0)
    )

    # Flatten column names
    snv_counts.columns = ['sweeps', 'total']
    snv_counts = snv_counts.reset_index()

    # Compute non_sweeps and ensure integer types
    snv_counts['non_sweeps'] = snv_counts['total'] - snv_counts['sweeps']
    snv_counts[['sweeps', 'total', 'non_sweeps']] = snv_counts[['sweeps', 'total', 'non_sweeps']].astype(int)

    # Order columns
    snv_counts = snv_counts[['patient_id', 'bin', 'group', 'sweeps', 'non_sweeps', 'total']]

    unique_pairs = snv_counts[['patient_id', 'bin']].drop_duplicates().shape[0]
    logger.info(f"SNV counts calculated for {unique_pairs} unique (patient_id, bin) combinations")

    return snv_counts

def create_frequency_histogram(snvs: pd.DataFrame, output_dir: Path) -> None:
    """Create frequency histogram using available frequency columns."""
    logger = logging.getLogger(__name__)
    logger.info("Creating frequency histogram...")
    
    # Use freq_range as the frequency measure (max_freq - min_freq)
    freq_data = snvs['freq_range'].dropna()
    
    if freq_data.empty:
        logger.warning("No frequency data available for histogram")
        return
    
    # Create histogram
    hist, bins = pd.cut(freq_data, bins=50, retbins=True)
    hist_counts = hist.value_counts().sort_index()
    
    # Create histogram DataFrame
    hist_df = pd.DataFrame({
        'bin_start': bins[:-1],
        'bin_end': bins[1:],
        'count': hist_counts.values,
        'bin_center': (bins[:-1] + bins[1:]) / 2
    })
    
    # Save histogram
    hist_path = output_dir / 'freq_range_histogram.csv'
    hist_df.to_csv(hist_path, index=False)
    logger.info(f"Frequency histogram saved to {hist_path}")


def create_frequency_histogram_by_group(snvs: pd.DataFrame, meta: pd.DataFrame, output_dir: Path) -> None:
    """Create group-wise frequency histogram in long format (group as a column)."""
    logger = logging.getLogger(__name__)
    logger.info("Creating frequency histogram by group (long format)...")

    # Join to get group
    snvs_with_group = snvs.merge(meta[['patient_id', 'group']], on='patient_id', how='left')

    if 'freq_range' not in snvs_with_group.columns:
        logger.warning("freq_range column missing; skipping group-wise histogram")
        return

    # Cut bins globally for consistent bin edges across groups
    freq_series = snvs_with_group['freq_range'].dropna()
    if freq_series.empty:
        logger.warning("No frequency data available for group-wise histogram")
        return

    _, bins = pd.cut(freq_series, bins=50, retbins=True)

    # Compute per-group histogram using the same bins
    rows = []
    for grp, df_g in snvs_with_group.groupby('group', dropna=True):
        if df_g.empty:
            continue
        cats = pd.cut(df_g['freq_range'].dropna(), bins=bins)
        counts = cats.value_counts().sort_index()
        # Build long-form rows
        for i in range(len(bins) - 1):
            bin_start = bins[i]
            bin_end = bins[i + 1]
            count = int(counts.iloc[i]) if i < len(counts) else 0
            rows.append({
                'group': grp,
                'bin_start': float(bin_start),
                'bin_end': float(bin_end),
                'bin_center': float((bin_start + bin_end) / 2),
                'count': count,
            })

    hist_long = pd.DataFrame(rows)
    out_path = output_dir / 'freq_range_histogram_by_group.csv'
    hist_long.to_csv(out_path, index=False)
    logger.info(f"Group-wise frequency histogram saved to {out_path}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Calculate SNV bin statistics')
    parser.add_argument('--input-dir', type=Path, required=True,
                       help='Directory containing all_snvs parquet files')
    parser.add_argument('--output-dir', type=Path, required=True,
                       help='Output directory for results')
    parser.add_argument('--meta-file', type=Path, required=True,
                       help='Metadata file path')
    
    args = parser.parse_args()
    
    # Setup
    logger = setup_logging()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load data
        meta = load_metadata(args.meta_file)
        snvs = load_snvs(args.input_dir)
        
        # Calculate SNVs per bin
        snv_counts = calculate_snvs_per_bin(snvs, meta)
        
        # Create frequency histogram (overall and by group)
        create_frequency_histogram(snvs, args.output_dir)
        create_frequency_histogram_by_group(snvs, meta, args.output_dir)
        
        # Save SNV counts
        snv_counts_path = args.output_dir / 'snv_counts_per_bin.csv'
        snv_counts.to_csv(snv_counts_path, index=False)
        logger.info(f"SNV counts per bin saved to {snv_counts_path}")
        
        # Save summary statistics
        summary = {
            'total_snvs': int(len(snvs)),
            'unique_patients': int(snvs['patient_id'].nunique()),
            'unique_bins': int(snvs['bin'].nunique()),
            'unique_patient_bin_pairs': int(snv_counts[['patient_id','bin']].drop_duplicates().shape[0]),
            'total_sweeps': int(snvs['is_sweep'].sum()),
            'sweep_percentage': float((snvs['is_sweep'].sum() / len(snvs)) * 100)
        }
        
        summary_path = args.output_dir / 'summary_stats.json'
        import json
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Summary statistics saved to {summary_path}")
        
        logger.info("All done!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
