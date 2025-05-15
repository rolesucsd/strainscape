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

# ────────── logging ──────────
def setup_logger(logfile=None):
    fmt = "%(asctime)s %(levelname)s  %(message)s"
    hdlr = logging.FileHandler(logfile) if logfile else logging.StreamHandler(sys.stderr)
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[hdlr])
    return logging.getLogger("trend")

log = setup_logger()


# ────────── main routine ──────────
def calc_trends_fast(muts_f  : Path,
                     meta_f  : Path,
                     out_f   : Path,
                     min_slope=0.01,
                     p_thr    =0.05):

    log.info("Load metadata")
    meta = (pd.read_csv(meta_f, usecols=["External.ID", "week_num"])
              .rename(columns={"External.ID":"Sample"})
              .drop_duplicates("Sample"))
    meta["week_num"] = meta["week_num"].astype(int)

    log.info("Load mutations")
    muts = pd.read_csv(muts_f, sep="\t",
                       usecols=["scaffold","position","Sample",
                                "ref_base","position_coverage",
                                "A","C","G","T"])

    # ── join + compute alt-allele frequencies
    muts = muts.merge(meta, on="Sample", how="inner")
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
        sub  = muts.loc[mask, ["scaffold","position","week_num"]].copy()
        sub["nucleotide"] = base
        sub["frequency"]  = freqs[mask,i]
        rows.append(sub)
    long = pd.concat(rows, ignore_index=True)

    # ── pivot to wide
    wide = (long.pivot_table(index=["scaffold","position","nucleotide"],
                             columns="week_num",
                             values="frequency",
                             aggfunc="first")
                 .sort_index(axis=1))
    week_cols = wide.columns.to_numpy(int)
    mat = wide.to_numpy(float)                       # (n_sites, n_weeks)

    # Require ≥3 valid time-points & some variance
    valid_mask = np.isfinite(mat).sum(1) >= 3
    mat = mat[valid_mask]
    wide = wide.iloc[valid_mask]
    var  = np.nanvar(mat, axis=1)
    mat  = mat[var>0]
    wide = wide.iloc[var>0]

    log.info(f"Sites for regression: {len(wide):,}")

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

    # ── filter by slope & p ──
    out = res[(np.abs(res["slope"])>=min_slope) & (res["p_value"]<=p_thr)]
    log.info(f"Significant: {len(out):,}")
    out.to_csv(out_f, sep="\t", index=False)
    log.info(f"Saved → {out_f}")


# ───────────── CLI ─────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mutation_file",  required=True)
    ap.add_argument("--metadata_file",  required=True)
    ap.add_argument("--output_file",    required=True)
    ap.add_argument("--min_slope",  type=float, default=0.01)
    ap.add_argument("--p_value",    type=float, default=0.05)
    ap.add_argument("--log_file")
    args = ap.parse_args()

    if args.log_file:
        log = setup_logger(args.log_file)

    calc_trends_fast(Path(args.mutation_file),
                     Path(args.metadata_file),
                     Path(args.output_file),
                     min_slope=args.min_slope,
                     p_thr=args.p_value)
