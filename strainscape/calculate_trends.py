#!/usr/bin/env python3
"""
Fast trend detection for SNV frequencies.

Strategy
────────
1.  Join mutations ⇢ metadata (adds `week_num`), drop rows without week.
2.  Build a *wide* matrix:  index = (scaffold, position, nucleotide),
    columns = week numbers, values = alt-allele frequency.
3.  Filter rows with < 3 time-points or zero variance.
4.  Use closed-form ordinary-least-squares to compute slope, r², p-value
    for every row at once (NumPy).
"""

from pathlib import Path
import logging, argparse, sys
import numpy as np, pandas as pd
from scipy.stats import t
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

# ────────── main routine ──────────
def calc_trends_fast(muts_f  : Path,
                     out_f   : Path,
                     min_slope=0.00,
                     p_thr    =1,
                     min_freq_range=0.2):
    """
    Identify nucleotide changes with significant frequency trends over time.
    
    PROCESS OVERVIEW:
    1. For each genomic position, calculate frequencies of A, C, G, T nucleotides
    2. Create separate rows for each alternative nucleotide (non-reference)
    3. Perform OLS regression to find which nucleotides show significant trends
    4. The 'new_base' in output represents the nucleotide with significant trend
    
    KEY INSIGHT: The "new_base" is determined by which nucleotide shows the most
    significant frequency change over time, not by comparing to a reference genome.
    This identifies evolutionary changes or selective pressures.
    """

    logger.info("Load mutations")
    muts = pd.read_csv(muts_f, sep="\t",
                       usecols=["scaffold","position","Sample",
                                "ref_base","position_coverage",
                                "week_num","A","C","G","T"])

    # STEP 1: Calculate allele frequencies for each nucleotide at each position
    bases = np.array(["A","C","G","T"])
    base_idx = muts["ref_base"].map({b:i for i,b in enumerate(bases)}).to_numpy()
    cov = muts["position_coverage"].to_numpy(float) # total coverage at position
    counts = muts[bases].to_numpy(float) # raw counts of each nucleotide (A, C, G, T)
    freqs  = counts / cov[:,None] # (n,4) # frequency of each nucleotide [0-1]

    # STEP 2: Create separate rows for each alternative nucleotide (non-reference)
    # This is the key step: for each position, we create rows for A, C, G, T
    # but only keep the ones that are different from the reference base
    # This allows us to test which alternative nucleotides show significant trends
    rows = []
    for i, base in enumerate(bases):
        # Only keep positions where this nucleotide is NOT the reference
        mask = base_idx != i  # True where ref_base != current base
        if not mask.any(): continue
        sub  = muts.loc[mask, ["scaffold","position","week_num", "ref_base", "position_coverage"]].copy()
        sub["nucleotide"] = base  # This becomes the "new_base" in output
        sub["frequency"]  = freqs[mask,i]  # Frequency of this nucleotide over time
        rows.append(sub)
    long = pd.concat(rows, ignore_index=True)

    # STEP 2.5: Write long format output with desired columns
    # Add position_coverage and calculate depth for each row
    long_with_coverage = []
    for _, row in long.iterrows():
        long_with_coverage.append({
            'bin_id': row['scaffold'].split('|')[0] if '|' in row['scaffold'] else row['scaffold'],
            'scaffold': row['scaffold'],
            'position': row['position'],
            'ref': row['ref_base'],
            'alt': row['nucleotide'],
            'af': row['frequency'],
            'coverage': row['position_coverage'],
            'time': row['week_num']
        })
    
    long_df = pd.DataFrame(long_with_coverage)
    long_out_f = out_f.parent / f"{out_f.stem}_long.tsv"
    long_df.to_csv(long_out_f, sep="\t", index=False)
    logger.info(f"Saved long format → {long_out_f}")

    # STEP 3: Pivot to wide format for regression analysis
    # Each row represents one (scaffold, position, nucleotide) combination
    # Columns are timepoints (weeks), values are allele frequencies
    wide = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                             columns="week_num",
                             values="frequency",
                             aggfunc="first")
                 .sort_index(axis=1))
    week_cols = wide.columns.to_numpy(int)  # Time points (weeks)
    mat = wide.to_numpy(float)              # (n_sites, n_weeks) frequency matrix

    # STEP 4: Filter for valid regression candidates
    # Require ≥3 valid time-points & some variance across time
    valid_mask = np.isfinite(mat).sum(1) >= 3  # At least 3 timepoints with data
    mat = mat[valid_mask]
    wide = wide.iloc[valid_mask]
    var  = np.nanvar(mat, axis=1)  # Variance across timepoints
    mat  = mat[var>0]              # Only keep sites with temporal variance
    wide = wide.iloc[var>0]

    logger.info(f"Sites for regression: {len(wide):,}")

    # STEP 5: Perform OLS regression for each nucleotide at each position
    # y = a + b*x where y=frequency, x=time, b=slope (trend)
    x   = week_cols.astype(float)  # Time points
    x_c = x - x.mean()             # Centered time points
    denom = (x_c**2).sum()         # Denominator for slope calculation
    b    = np.nansum(mat * x_c, axis=1) / denom  # SLOPE: rate of frequency change
    a    = np.nanmean(mat, axis=1)               # Intercept (mean frequency)
    y_hat = a[:,None] + b[:,None]*x              # Predicted frequencies
    ss_tot = np.nansum((mat - mat.mean(1,keepdims=True))**2, axis=1)  # Total sum of squares
    ss_res = np.nansum((mat - y_hat)**2, axis=1)                      # Residual sum of squares
    r2   = 1 - ss_res/ss_tot  # R-squared: proportion of variance explained

    # STEP 6: Calculate statistical significance of trends
    n    = np.isfinite(mat).sum(1)  # Number of valid timepoints per site
    se   = np.sqrt(ss_res / (n-2)) / np.sqrt(denom)  # Standard error of slope
    tval = b / se                   # t-statistic for slope
    pval = 2 * t.sf(np.abs(tval), df=n-2)  # Two-sided p-value

    # STEP 7: Calculate summary statistics for each nucleotide trend
    minf = np.nanmin(mat, axis=1)      # Minimum frequency across time
    maxf = np.nanmax(mat, axis=1)      # Maximum frequency across time

    # STEP 8: Assemble final results
    res = wide.reset_index()
    res["slope"]      = b              # Rate of frequency change (positive = increasing, negative = decreasing)
    res["p_value"]    = pval           # Statistical significance of trend
    res["r_squared"]  = r2             # How well the linear trend fits the data
    res["min_freq"]   = minf           # Lowest frequency observed
    res["max_freq"]   = maxf           # Highest frequency observed
    res["freq_range"] = maxf - minf    # Total range of frequencies
    res["mean_freq"]  = np.nanmean(mat, axis=1)  # Average frequency across time
    res = res.rename(columns={"scaffold":"scaffold",
                              "position":"position",
                              "nucleotide":"new_base"})  # This is the key: nucleotide becomes new_base

    # STEP 9: Filter for significant trends
    # Keep only trends that meet our criteria:
    # - Absolute slope >= min_slope (minimum rate of change)
    # - p-value <= p_thr (statistical significance)
    # - Frequency range >= min_freq_range (minimum variation over time)
    out = res[(np.abs(res["slope"])>=min_slope) & (res["p_value"]<=p_thr) & (res["freq_range"]>=min_freq_range)]
    logger.info(f"Significant: {len(out):,}")
    
    out.to_csv(out_f, sep="\t", index=False)
    logger.info(f"Saved → {out_f}")


# ───────────── CLI ─────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mutation_file",  required=True)
    ap.add_argument("--output_file",    required=True)
    ap.add_argument("--min_slope",  type=float, default=0.0)
    ap.add_argument("--p_value",    type=float, default=1.0)
    ap.add_argument("--min_freq_range",     type=float, default=0.2)
    ap.add_argument("--log_file", required=False)
    args = ap.parse_args()

    if args.log_file:
        setup_logging(args.log_file)

    calc_trends_fast(Path(args.mutation_file),
                     Path(args.output_file),
                     min_slope=args.min_slope,
                     p_thr=args.p_value,
                     min_freq_range=args.min_freq_range)
