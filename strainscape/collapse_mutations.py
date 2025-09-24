#!/usr/bin/env python3
"""
fast_collapse_mutations_birch.py  (≤200 lines)

Per-bin clustering of mutation trajectories with BIRCH (Euclidean).
- Filters rows by dynamic range (default ≥0.25) and min observed points.
- Mean-imputes remaining NAs per timepoint (within bin) for distance calcs.
- BIRCH clusters with radius `--threshold` (Euclidean).
- Optional final global reclustering of BIRCH subclusters.
- Parallel across bins; writes per-bin and combined summaries.
"""

import argparse, os, re
import logging
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Any
import numpy as np
import pandas as pd
from sklearn.cluster import Birch
from sklearn.cluster import KMeans

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ---------- helpers ----------
def detect_traj_cols(df: pd.DataFrame) -> List[str]:
    """Timepoint columns = all columns that look numeric or hold numeric values."""
    logger.info(f"Detecting trajectory columns from {len(df.columns)} total columns")
    
    out = []
    regex_matches = 0
    
    for c in df.columns:
        if re.fullmatch(r'[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?', str(c)):
            out.append(c)
            regex_matches += 1
            continue
    
    logger.info(f"Detected {len(out)} trajectory columns: {regex_matches} by regex")
    logger.info(f"First 10 trajectory columns: {out[:10]}")
    return out


# ---------- core per-bin ----------
def cluster_bin_birch(
    df_bin: pd.DataFrame,
    traj_cols: List[str],
    freq_range_thresh: float,
    min_pts: int,
    threshold: float,
    branching_factor: int,
    random_state: int
) -> pd.DataFrame:
    logger.info(f"Processing bin with {len(df_bin)} mutations and {len(traj_cols)} timepoints")
    
    X = df_bin[traj_cols].to_numpy(dtype=np.float32)
    mask = np.isfinite(X)
    logger.info(f"Data shape: {X.shape}, finite values: {mask.sum()}/{mask.size} ({mask.sum()/mask.size*100:.1f}%)")

    # pre-filter: require enough observations & dynamic range on observed points
    obs = np.where(mask, X, np.nan)
    row_min = np.nanmin(obs, axis=1)
    row_max = np.nanmax(obs, axis=1)
    n_obs   = mask.sum(axis=1)
        
    keep = (n_obs >= min_pts) & np.isfinite(row_min) & np.isfinite(row_max) & ((row_max - row_min) >= freq_range_thresh)
    logger.info(f"After filtering: {keep.sum()}/{len(keep)} mutations kept ({keep.sum()/len(keep)*100:.1f}%)")
    
    if not np.any(keep): 
        logger.warning("No mutations passed filtering criteria")
        return df_bin.iloc[0:0].assign(cluster=[])

    X = X[keep]; mask = mask[keep]
    kept_idx = np.where(keep)[0]

    # mean-impute per timepoint within bin (required for distances)
    col_means = np.nanmean(np.where(mask, X, np.nan), axis=0)
    # guard: if a column is all-NaN, set mean to 0
    col_means = np.where(np.isfinite(col_means), col_means, 0.0).astype(np.float32)
    X_imp = np.where(mask, X, col_means)
    
    # BIRCH pass
    logger.info(f"Running BIRCH clustering with threshold={threshold}, branching_factor={branching_factor}")
    thr = (len(traj_cols) * threshold)
    birch = Birch(threshold=thr, branching_factor=branching_factor, n_clusters=None)
    birch.fit(X_imp)
    labels = birch.predict(X_imp)
    
    n_clusters = len(np.unique(labels))
    logger.info(f"BIRCH created {n_clusters} clusters")
    
    # Log cluster size distribution
    unique, counts = np.unique(labels, return_counts=True)
    logger.info(f"Cluster sizes - min: {counts.min()}, max: {counts.max()}, mean: {counts.mean():.1f}")

    out = df_bin.iloc[kept_idx].copy()
    out["cluster"] = labels.astype(int)
    return out

def process_one_bin(bin_id: Any, df_bin: pd.DataFrame, traj_cols: List[str], args) -> pd.DataFrame:
    logger.info(f"=== Processing bin: {bin_id} ===")
    clustered = cluster_bin_birch(
        df_bin, traj_cols,
        freq_range_thresh=args.freq_range,
        min_pts=args.min_pts,
        threshold=args.threshold,
        branching_factor=args.branching_factor,
        random_state=args.seed
    )
    if clustered.empty: 
        logger.warning(f"Bin {bin_id} produced no clusters")
        return clustered
    
    logger.info(f"Summarizing clusters for bin {bin_id}")
    output_file = os.path.join(args.outdir, f"{bin_id}_clusters.tsv")
    clustered.to_csv(output_file, sep="\t", index=False)
    return clustered

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Fast per-bin BIRCH clustering of mutation trajectories (Euclidean).")
    ap.add_argument("input", help="Wide TSV/CSV from analyze_mutations.py")
    ap.add_argument("outdir", help="Output directory")
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--freq-range", type=float, default=0.25, help="Min allele-frequency range to keep a row")
    ap.add_argument("--min-pts", type=int, default=3, help="Min non-NA timepoints required per row")
    ap.add_argument("--threshold", type=float, default=0.05, help="BIRCH radius (Euclidean) for cluster formation")
    ap.add_argument("--branching-factor", type=int, default=50)
    ap.add_argument("--final-clusters", type=int, default=0, help="Optional k for final clustering of BIRCH CFs (0=off)")
    ap.add_argument("--workers", type=int, default=max(1, cpu_count()//2))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    logger.info("=== Starting BIRCH Clustering ===")
    logger.info(f"Parameters: freq_range={args.freq_range}, min_pts={args.min_pts}, threshold={args.threshold}")
    logger.info(f"Branching factor: {args.branching_factor}, Workers: {args.workers}")

    os.makedirs(args.outdir, exist_ok=True)
    logger.info(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input, sep=args.sep, low_memory=False)
    logger.info(f"Loaded {len(df):,} rows and {len(df.columns)} columns")

    traj_cols = detect_traj_cols(df)
    if not traj_cols: 
        logger.error("No timepoint columns detected!")
        raise SystemExit("No timepoint columns detected.")
    
    logger.info("Converting trajectory columns to numeric...")
    for c in traj_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float32")

    groups = list(df.groupby("bin", sort=False))
    logger.info(f"Found {len(groups)} bins to process")
    logger.info(f"Bin names: {[g[0] for g in groups[:5]]}{'...' if len(groups) > 5 else ''}")
    
    worker = partial(process_one_bin, traj_cols=traj_cols, args=args)
    logger.info(f"Starting processing with {args.workers} workers")
    
    results = (Pool(args.workers).starmap(worker, groups)
               if args.workers > 1 and len(groups) > 1
               else [worker(b, g) for b, g in groups])

    summaries = [r for r in results if r is not None and not r.empty]
    logger.info(f"Processing complete. {len(summaries)} bins produced results")
    
    if summaries:
        combined_file = os.path.join(args.outdir, "all_bins_clusters.tsv")
        combined_df = pd.concat(summaries, ignore_index=True)
        combined_df.to_csv(combined_file, sep="\t", index=False)
        logger.info(f"Saved combined results: {len(combined_df):,} cluster summaries to {combined_file}")
    else:
        logger.warning("No clusters were produced!")

if __name__ == "__main__":
    main()
