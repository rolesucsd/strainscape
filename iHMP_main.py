#!/usr/bin/env python3
"""
This script processes sample directories using pre-loaded metadata for SNP analysis in the iHMP dataset.
It applies filtering, trend analysis, clustering, mapping, and mutation typing in a modular pipeline.
"""

import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
from typing import Optional
import os
from statsmodels.stats.multitest import multipletests
import glob
import logging
from tqdm import tqdm

# Import custom functions from modules (assumed to exist)
from iHMP_trends import calculate_trend
from iHMP_map import read_annotation2, batch_map_genes
from iHMP_mutation import get_mutation_type, visualize_genic_vs_intergenic_and_mutations
from iHMP_filter import read_input, filter_snv, normalize_counts_vec
from iHMP_cluster import add_selection_coeff, cluster_snps_combined
from iHMP_trajectory import prepare_trajectory

# ------------------------------------------------------------------------------
# Logging Setup
# ------------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Log Python version for debugging
logger.info(f"Running with Python version: {sys.version}")

# ------------------------------------------------------------------------------
# Global Constants & Metadata Loading
# ------------------------------------------------------------------------------
METADATA_FILE = "metadata/hmp2_metadata_2018-08-20.csv"
GENES_FILE = "../reference/bakta/wolr2_reference.tsv"
SPECIES_FILE = "../reference/nucleotide_div_filter.txt"
FASTA_FILE = "../reference/wolr2_reference.fna"
TREES_DIR = "../reference/trees"
BASE_PATTERN = "../instrain_output/*/"

# Load metadata files once globally
try:
    metadata1 = pd.read_csv("../reference/wolr2_reference.stb", sep='\t')
    metadata2 = pd.read_csv("../reference/all_unique_genomes_assembly.txt", sep='\t', usecols=[0, 1])
    metadata3 = pd.read_csv("../reference/lineages.txt", sep='\t')
    FINAL_MERGED = pd.merge(pd.merge(metadata1, metadata2, on="Assembly", how="inner"),
                            metadata3, on="Genome", how="inner")
    SEQUENCES = {record.id: str(record.seq) for record in SeqIO.parse(FASTA_FILE, "fasta")}
except FileNotFoundError as e:
    logger.error(f"Failed to load metadata or FASTA file: {e}")
    raise

# ------------------------------------------------------------------------------
# Helper Function for I/O Simplification
# ------------------------------------------------------------------------------
def load_or_process(output_file: str, process_func, *args) -> Optional[pd.DataFrame]:
    """
    Load an existing file or process and save data if it doesn't exist.

    Args:
        output_file (str): Path to the output file.
        process_func (callable): Function to generate data if file is missing.
        *args: Arguments to pass to process_func.

    Returns:
        Optional[pd.DataFrame]: Processed or loaded DataFrame, or None if processing fails.
    """
    if not os.path.exists(output_file):
        logger.info(f"Creating {output_file}")
        result = process_func(*args)
        if result is not None and not result.empty:
            result.to_csv(output_file, sep='\t', index=False)
            return result
        logger.warning(f"No data generated for {output_file}")
        return None
    logger.info(f"Reading {output_file}")
    return pd.read_csv(output_file, sep='\t', low_memory=False)

# ------------------------------------------------------------------------------
# Process a Single Sample Directory
# ------------------------------------------------------------------------------
def process_sample_directory(file_dir: str, merged_metadata: pd.DataFrame) -> None:
    """
    Process a single sample directory through the SNP analysis pipeline.

    Args:
        file_dir (str): Directory path for the sample (must end with a slash).
        merged_metadata (pd.DataFrame): Pre-loaded merged metadata DataFrame.

    Returns:
        None: Writes output files to disk or skips if they exist.
    """
    # Define input and output file paths
    snv_master_file = os.path.join(file_dir, "combined_SNV_info.tsv")
    output_files = {
        1: os.path.join(file_dir, "new_SNV_filtered.txt"),
        2: os.path.join(file_dir, "new_SNV_filtered_trend.txt"),
        3: os.path.join(file_dir, "new_SNV_filtered_trend_trajectory.txt"),
        4: os.path.join(file_dir, "new_SNV_filtered_trend_trajectory_mapped.txt"),
        5: os.path.join(file_dir, "new_SNV_filtered_trend_trajectory_mapped_mutation.txt")
    }

    # Input validation
    if not os.path.exists(snv_master_file):
        logger.error(f"SNV master file not found: {snv_master_file}")
        return
    if os.path.getsize(snv_master_file) == 0:
        logger.warning(f"SNV master file is empty: {snv_master_file}")
        return

    # Step 1: Filter and normalize SNVs
    SNV_normalized = load_or_process(
        output_files[1],
        lambda: normalize_counts_vec(filter_snv(read_input(snv_master_file, SPECIES_FILE, METADATA_FILE), save_fig=file_dir))
    )
    if SNV_normalized is None:
        return
    logger.info(SNV_normalized.head())


    # Step 2: Calculate trends
    SNV_trend = load_or_process(
        output_files[2],
        lambda: SNV_normalized.merge(calculate_trend(SNV_normalized, n_jobs=4), on=['scaffold', 'position'], how='right')
    )
    if SNV_trend is None or SNV_trend.empty:
        logger.info(f"No SNPs after quality filtering in {file_dir}")
        return

    # Step 3: Prepare trajectory
    SNV_trajectory = load_or_process(
        output_files[3],
        lambda: pd.merge(prepare_trajectory(SNV_trend), merged_metadata, on='Chromosome', how="left")
    )
    if SNV_trajectory is None:
        return

    # Step 4: No filtering!
    if not os.path.exists(output_files[4]):
        logger.info(f"Creating {output_files[4]}")
        if 'OLS_pvalue' not in SNV_trajectory.columns:
            logger.error(f"OLS_pvalue column missing in {file_dir}")
            return
        mask_valid = SNV_trajectory['OLS_pvalue'].notnull()
        if mask_valid.sum() > 0:
            _, adj_pvals, _, _ = multipletests(SNV_trajectory.loc[mask_valid, 'OLS_pvalue'], method='fdr_bh')
            SNV_trajectory.loc[mask_valid, 'OLS_pvalue_adjusted'] = adj_pvals
        else:
            SNV_trajectory['OLS_pvalue_adjusted'] = np.nan
            logger.warning(f"No valid p-values in {file_dir}")

    # Filter out all rows with OLS_pvalue >= 0.1
    before_n = len(SNV_trajectory)
    SNV_trajectory = SNV_trajectory[SNV_trajectory['OLS_pvalue'] < 0.1].copy()
    after_n  = len(SNV_trajectory)
    logger.info(f"Filtered out {before_n - after_n} rows with OLS_pvalue >= 0.1")

    # Step 4 cont: Map genes
    SNV_mapped = load_or_process(
        output_files[4],
        lambda: process_mapping(SNV_trajectory, GENES_FILE, TREES_DIR)
    )
    if SNV_mapped is None:
        return

    # Step 5: Add mutation types and visualize
    if not os.path.exists(output_files[5]):
        logger.info(f"Creating {output_files[5]}")
        SNV_mapped['Mutation_Type'] = SNV_mapped.apply(get_mutation_type, axis=1, sequences=SEQUENCES)
        visualize_genic_vs_intergenic_and_mutations(SNV_mapped, save_fig=file_dir)
        SNV_mapped.to_csv(output_files[5], sep='\t', index=False)
    else:
        logger.info(f"Processing finished in {file_dir}")

def process_mapping(SNV_trajectory: pd.DataFrame, genes_file: str, trees_dir: str) -> pd.DataFrame:
    """
    Helper function to filter and map genes to clustered SNVs.

    Args:
        SNV_trajectory (pd.DataFrame): Clustered SNV data.
        genes_file (str): Path to gene annotation file.
        trees_dir (str): Directory containing interval trees.

    Returns:
        pd.DataFrame: Mapped SNV data.
    """
    genes = read_annotation2(genes_file)
    snv_mapped = batch_map_genes(SNV_trajectory, genes, trees_dir)
    return pd.merge(SNV_trajectory, snv_mapped, on=['Chromosome', 'Position'], how='left')

# ------------------------------------------------------------------------------
# Main Function
# ------------------------------------------------------------------------------
def main() -> None:
    """Main function to process all sample directories with a progress bar."""
    subdirs = sorted(glob.glob(BASE_PATTERN))
    for file_dir in tqdm(subdirs, desc="Processing directories"):
        file_dir = os.path.abspath(file_dir) + os.sep
        output6 = os.path.join(file_dir, "SNV_filtered_trend_cluster_mapped_mutation.txt")
        logger.info(f"\nProcessing sample directory: {file_dir}")
        if not os.path.exists(output6):
            try:
                process_sample_directory(file_dir, FINAL_MERGED)
            except Exception as e:
                logger.error(f"Error processing {file_dir}: {e}")
        else:
            logger.info(f"Processing finished for {file_dir}")

if __name__ == "__main__":
    main()