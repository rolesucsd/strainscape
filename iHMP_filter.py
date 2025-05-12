import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.metrics import r2_score
from joblib import Parallel, delayed
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def read_input(snv_master_file, species_master_file, metadata_file):
    """
    Reads input SNV master and metadata files, then merges them on 'External.ID'.

    Parameters:
        snv_master_file (str): Path to the SNV master file (TSV format).
        metadata_file (str): Path to the metadata file (CSV format).

    Returns:
        pd.DataFrame: Merged dataframe containing SNV and metadata information.
    """
    SNV_master = pd.read_csv(snv_master_file, low_memory=False, sep="\t")
    metadata = pd.read_csv(metadata_file, low_memory=False)
    
    # Subset metadata to include only relevant columns
    metadata = metadata[['External.ID', 'week_num', 'Participant ID']].drop_duplicates()
    
    # Merge SNV master with metadata
    SNV_wide = SNV_master.merge(metadata, on='External.ID', how='left')

    species_df = pd.read_csv(species_master_file, sep='\t')

    # Build a small dataframe of valid (scaffold, Sample) pairs
    valid_pairs = species_df[['scaffold','Participant ID']].drop_duplicates()

    # Inner‐merge will keep only rows in SNV_wide that match those pairs
    SNV_wide_filtered = SNV_wide.merge(
        valid_pairs,
        on=['scaffold','Participant ID'],
        how='inner'
        )
    
    return SNV_wide_filtered

def plot_coverage_distribution(coverage_col, coverage_threshold, steps, counts, save_fig=None):
    """
    Plots the coverage distribution for SNPs and saves the data used for plotting.

    Parameters:
        coverage_col (pd.Series): Column containing coverage values.
        coverage_threshold (int): Threshold for filtering SNPs based on coverage.
        steps (list of str): Names of filtering steps.
        counts (list of int): Counts of SNPs after each filtering step.
        save_fig (str): Directory to save the plot and data file. If None, the plot is displayed.
    """
    # Plot coverage distribution BEFORE filtering
    fig = plt.figure(figsize=(10, 6))
    
    # Coverage histogram
    gs = GridSpec(1, 2, width_ratios=[7, 3])

    ax1 = fig.add_subplot(gs[0])

    # Log-spaced bins from 1 to 4000
    bins = np.logspace(np.log10(1), np.log10(4000), 20)

    sns.histplot(coverage_col, bins=bins, color='steelblue', ax=ax1)
    ax1.set_xscale('log')

    ax1.axvline(coverage_threshold, color='red', linestyle='--', label=f"Coverage Threshold = {coverage_threshold}")
    ax1.legend()
    ax1.set_title("Coverage Distribution (Before Filtering)")
    ax1.set_xlabel("Coverage")
    ax1.set_ylabel("Number of SNPs")

    ax2 = fig.add_subplot(gs[1])

    sns.barplot(x=steps, y=counts, palette="viridis", ax=ax2)
    ax2.set_title("SNP Counts After Each Filter")
    for i, c in enumerate(counts):
        ax2.text(i, c + 0.01 * c, str(c), ha='center', va='bottom')
    for label in ax2.get_xticklabels():
        label.set_rotation(45)
    ax2.set_ylabel("Number of SNPs")

    fig.tight_layout()

    if save_fig:
        # Save the figure
        fig.savefig(os.path.join(save_fig, "snp_filtering_steps.png"), dpi=300)

        # Save the data used for plotting
        data = pd.DataFrame({
            "Step": steps,
            "SNP_Count": counts
        })
        data.to_csv(os.path.join(save_fig, "snp_filtering_data.csv"), index=False)
    else:
        plt.show()

def filter_snv(SNV_df, coverage_col='position_coverage', coverage_threshold=10, 
               groupby_cols=('scaffold', 'position'), group_min_size=5, 
               save_fig=None):
    """
    Filters SNVs based on coverage, group size, and (optionally) allele frequency values.
    
    In addition to filtering out SNPs with coverage below a threshold and groups
    with fewer than a minimum number of records, this function removes any rows
    for which the values in the specified frequency columns (after ignoring NaNs)
    are all identical and equal to 0 or all identical and equal to 100.
    
    Parameters:
        SNV_df (pd.DataFrame): DataFrame containing SNP data.
        coverage_col (str): Column name for coverage values.
        coverage_threshold (int): Minimum coverage to retain SNPs.
        groupby_cols (tuple of str): Columns to group by for filtering.
        group_min_size (int): Minimum group size to retain SNPs.
        save_fig (str, optional): Path to save figures. Defaults to None.
        
    Returns:
        pd.DataFrame: Filtered and normalized SNP DataFrame.
    """
    total_snps_before = len(SNV_df)
    
    # Filter by coverage threshold.
    SNV_df = SNV_df[SNV_df[coverage_col] >= coverage_threshold]
    total_after_coverage = len(SNV_df)
    
    # Filter by group size.
    group_sizes = SNV_df.groupby(list(groupby_cols)).size()
    valid_positions = group_sizes[group_sizes >= group_min_size].index
    SNV_df = SNV_df.set_index(list(groupby_cols)).loc[valid_positions].reset_index()
    total_after_group = len(SNV_df)
    
    # Plot filtering steps (assuming plot_coverage_distribution is defined elsewhere).
    counts = [total_snps_before, total_after_coverage, total_after_group]
    steps = [f"Original", f"Coverage >= {coverage_threshold}", f"Group size >= {group_min_size}"]
    plot_coverage_distribution(SNV_df[coverage_col], coverage_threshold, steps, counts, save_fig)
    
    return SNV_df

def normalize_counts_vec(df):
    """
    Normalizes nucleotide counts (A, T, C, G) as percentages.

    Parameters:
        df (pd.DataFrame): Dataframe with nucleotide count columns (A, T, C, G).

    Returns:
        pd.DataFrame: Dataframe with normalized nucleotide percentages.
    """
    total = df[['A', 'T', 'C', 'G']].sum(axis=1)
    mask = total > 0
    df.loc[mask, ['A', 'T', 'C', 'G']] = df.loc[mask, ['A', 'T', 'C', 'G']].div(total[mask], axis=0)
    df.loc[~mask, ['A', 'T', 'C', 'G']] = 0
    return df

