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
                     min_slope=0.01,
                     p_thr    =0.05):

    logger.info("Load mutations")
    muts = pd.read_csv(muts_f, sep="\t",
                       usecols=["scaffold","position","Sample",
                                "ref_base","position_coverage",
                                "week_num","A","C","G","T"])

    bases = np.array(["A","C","G","T"])
    base_idx = muts["ref_base"].map({b:i for i,b in enumerate(bases)}).to_numpy()
    cov = muts["position_coverage"].to_numpy(float)
    counts = muts[bases].to_numpy(float)
    freqs  = counts / cov[:,None]                    # (n,4)

    # Melt-free: create one row per alt base ≠ ref
    rows = []
    for i, base in enumerate(bases):
        mask = base_idx != i
        if not mask.any(): continue
        sub  = muts.loc[mask, ["scaffold","position","week_num", "ref_base", "position_coverage"]].copy()
        sub["nucleotide"] = base
        sub["frequency"]  = freqs[mask,i]
        rows.append(sub)
    long = pd.concat(rows, ignore_index=True)

    # ── pivot to wide for frequencies
    wide = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                             columns="week_num",
                             values="frequency",
                             aggfunc="first")
                 .sort_index(axis=1))
    
    # ── pivot to wide for depth (aggregate: mean, std, min, max)
    depth_wide_mean = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                                       columns="week_num",
                                       values="position_coverage",
                                       aggfunc="mean")
                          .sort_index(axis=1))
    depth_wide_std = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                                      columns="week_num",
                                      values="position_coverage",
                                      aggfunc=lambda x: x.std() if len(x) > 1 else 0.0)
                         .sort_index(axis=1))
    depth_wide_min = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                                       columns="week_num",
                                       values="position_coverage",
                                       aggfunc="min")
                          .sort_index(axis=1))
    depth_wide_max = (long.pivot_table(index=["scaffold","position","nucleotide","ref_base"],
                                       columns="week_num",
                                       values="position_coverage",
                                       aggfunc="max")
                          .sort_index(axis=1))
    week_cols = wide.columns.to_numpy(int)
    mat = wide.to_numpy(float)                       # (n_sites, n_weeks)

    # Require ≥3 valid time-points & some variance
    valid_mask = np.isfinite(mat).sum(1) >= 3
    mat = mat[valid_mask]
    wide = wide.iloc[valid_mask]
    # Apply same filtering to depth matrices
    depth_wide_mean = depth_wide_mean.iloc[valid_mask]
    depth_wide_std = depth_wide_std.iloc[valid_mask]
    depth_wide_min = depth_wide_min.iloc[valid_mask]
    depth_wide_max = depth_wide_max.iloc[valid_mask]
    
    var  = np.nanvar(mat, axis=1)
    mat  = mat[var>0]
    wide = wide.iloc[var>0]
    # Apply same filtering to depth matrices
    depth_wide_mean = depth_wide_mean.iloc[var>0]
    depth_wide_std = depth_wide_std.iloc[var>0]
    depth_wide_min = depth_wide_min.iloc[var>0]
    depth_wide_max = depth_wide_max.iloc[var>0]

    logger.info(f"Sites for regression: {len(wide):,}")

    # ── closed-form OLS y = a + b x ──
    x   = week_cols.astype(float)
    x_c = x - x.mean()
    denom = (x_c**2).sum()
    b    = np.nansum(mat * x_c, axis=1) / denom          # slope
    a    = np.nanmean(mat, axis=1)                       # intercept (unused)
    y_hat = a[:,None] + b[:,None]*x
    ss_tot = np.nansum((mat - mat.mean(1,keepdims=True))**2, axis=1)
    ss_res = np.nansum((mat - y_hat)**2, axis=1)
    r2   = 1 - ss_res/ss_tot

    # two-sided p-value for slope
    n    = np.isfinite(mat).sum(1)
    se   = np.sqrt(ss_res / (n-2)) / np.sqrt(denom)
    tval = b / se
    pval = 2 * t.sf(np.abs(tval), df=n-2)

    # min / max / range / mean
    minf = np.nanmin(mat, axis=1)
    maxf = np.nanmax(mat, axis=1)

    # ── assemble result
    res = wide.reset_index()
    res["slope"]      = b
    res["p_value"]    = pval
    res["r_squared"]  = r2
    res["min_freq"]   = minf
    res["max_freq"]   = maxf
    res["freq_range"] = maxf - minf
    res["mean_freq"]  = np.nanmean(mat, axis=1)
    res = res.rename(columns={"scaffold":"scaffold",
                              "position":"position",
                              "nucleotide":"new_base"})
    
    # Add depth statistics (aggregated across time points)
    depth_mat_mean = depth_wide_mean.to_numpy(float)
    depth_mat_std = depth_wide_std.to_numpy(float)
    depth_mat_min = depth_wide_min.to_numpy(float)
    depth_mat_max = depth_wide_max.to_numpy(float)
    
    res["depth_mean"] = np.nanmean(depth_mat_mean, axis=1)
    res["depth_std"] = np.nanstd(depth_mat_mean, axis=1)  # std of means across time
    res["depth_min"] = np.nanmin(depth_mat_min, axis=1)    # min of mins
    res["depth_max"] = np.nanmax(depth_mat_max, axis=1)    # max of maxs

    # ── filter by slope & p ──
    out = res[(np.abs(res["slope"])>=min_slope) & (res["p_value"]<=p_thr)]
    logger.info(f"Significant: {len(out):,}")
    out.to_csv(out_f, sep="\t", index=False)
    logger.info(f"Saved → {out_f}")


# ───────────── CLI ─────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mutation_file",  required=True)
    ap.add_argument("--output_file",    required=True)
    ap.add_argument("--min_slope",  type=float, default=0.01)
    ap.add_argument("--p_value",    type=float, default=0.05)
    ap.add_argument("--log_file", required=False)
    args = ap.parse_args()

    if args.log_file:
        setup_logging(args.log_file)

    calc_trends_fast(Path(args.mutation_file),
                     Path(args.output_file),
                     min_slope=args.min_slope,
                     p_thr=args.p_value)
