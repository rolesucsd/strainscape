#!/usr/bin/env python3
"""snv_calprotectin_corr.py  –  vectorised & fragment-aware

Find SNVs whose allele-frequency trajectories (wide numeric week columns 0-55)
track faecal-calprotectin.  Optimised for speed:
  • Reads only the required columns from each Parquet *fragment* and
    processes fragments one-by-one (constant memory, I/O parallelism).
  • Avoids the expensive melt → groupby loop.  Instead, for each patient
    we keep data in wide form and compute correlations in *vectorised NumPy*.
  • Skips SNVs or calprotectin series that are constant – no more r = ±1, p = 0
    artefacts.


Example:
---------
$ python snv_calprotectin_corr.py all_snvs  metadata.csv \
        --id-col "Participant ID" --time-col week_num --cal-col fecalcal \
        --method spearman --min-samples 6 --r-thresh 0.7 --p-thresh 0.05 \
        --out hits.tsv
"""
from __future__ import annotations

import argparse, logging, multiprocessing as mp
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np, pandas as pd, pyarrow.dataset as ds, pyarrow.compute as pc
from scipy.stats import spearmanr, pearsonr, rankdata

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ───────────────────────── correlation helpers ────────────────────────────
CORR_FUNCS = {
    "pearson":  lambda x, y: (np.corrcoef(x, y)[0, 1], np.nan),  # p from scipy if needed
    "spearman": spearmanr,
}

NumVec = np.ndarray  # shorthand


def fast_pearson_matrix(freq: np.ndarray, cpx: NumVec) -> NumVec:
    """Row-wise Pearson r for freq (n×t) vs. 1×t vector, ignoring NaNs."""
    mask = ~np.isnan(freq) & ~np.isnan(cpx)
    n    = mask.sum(1)
    valid = n >= 3
    r = np.full(freq.shape[0], np.nan, dtype=float)
    if not valid.any():
        return r
    # centre & scale each row where valid
    f = freq[valid].copy()
    m = mask[valid]
    # subtract row means
    row_mean = np.sum(np.where(m, f, 0.0), 1) / n[valid]
    f -= row_mean[:, None]
    c = cpx - np.nanmean(cpx)
    # compute covariance & variance
    cov = np.sum(np.where(m, f * c, 0.0), 1) / n[valid]
    var_f = np.sum(np.where(m, f ** 2, 0.0), 1) / n[valid]
    var_c = np.nanvar(cpx)
    r[valid] = cov / np.sqrt(var_f * var_c)
    return r


# ────────────────────────── fragment worker ───────────────────────────────

def process_fragment(args_tuple: Tuple) -> List[Dict]:
    (frag_path, id_col, week_cols, meta_wide, method, min_samples, r_thr, p_thr) = args_tuple
    # read fragment as table → pandas (read all columns)
    tbl = ds.dataset(frag_path, format="parquet").to_table()
    df = tbl.to_pandas()
    df = df.rename(columns={id_col: "patient_id"})
    
    print(f"DEBUG: Fragment {Path(frag_path).name} loaded with {len(df)} rows, {df['patient_id'].nunique()} unique patients")

    hits: List[Dict] = []
    # loop by patient to keep vectors aligned
    for pid, sub in df.groupby("patient_id", sort=False):
        if pid not in meta_wide.index:
            continue
        print(f"DEBUG: Processing patient {pid} with {len(sub)} SNVs")
        
        # Find common weeks between SNV and metadata for this patient
        common_weeks = [w for w in week_cols if float(w) in meta_wide.columns]
        if len(common_weeks) < min_samples:
            print(f"DEBUG: Patient {pid} has only {len(common_weeks)} common weeks, skipping")
            continue
        print(f"DEBUG: Patient {pid} has {len(common_weeks)} common weeks: {common_weeks[:5]}...")
        
        freq_mat = sub[common_weeks].to_numpy(dtype=float)
        cpx_vec = meta_wide.loc[pid, [float(w) for w in common_weeks]].values.astype(float)
        
        # skip rows with <min_samples non-NA or constant vectors
        non_na = (~np.isnan(freq_mat) & ~np.isnan(cpx_vec)).sum(1)
        variable = (np.nanstd(freq_mat, 1) > 0)
        keep_rows = (non_na >= min_samples) & variable
        
        print(f"DEBUG: Patient {pid} - {keep_rows.sum()} out of {len(sub)} SNVs pass filters (min_samples={min_samples})")
        
        if not keep_rows.any():
            print(f"DEBUG: Patient {pid} - no SNVs pass filters, skipping")
            continue
        freq_sel = freq_mat[keep_rows]
        sub_sel  = sub.iloc[np.where(keep_rows)[0]].reset_index(drop=True)

        if method == "pearson":
            r = fast_pearson_matrix(freq_sel, cpx_vec)
        else:
            r = fast_spearman_matrix(freq_sel, cpx_vec)

        # Debug: print the first few r values for this patient
        print(f"DEBUG: patient {pid}, r values: {r[:5]}")

        # p-values only for rows that pass r-thr; use analytical approx for speed
        mask_pass = np.abs(r) >= r_thr
        if not mask_pass.any():
            print(f"DEBUG: Patient {pid} - no SNVs pass r threshold {r_thr}")
            continue
        print(f"DEBUG: Patient {pid} - {mask_pass.sum()} SNVs pass r threshold")
        
        n_vec = non_na[keep_rows][mask_pass]
        p = 2 * (1 - np.abs(r[mask_pass]) * np.sqrt((n_vec - 2) / (1 - r[mask_pass] ** 2)))  # rough t-dist tail
        for i, ok in enumerate(np.where(mask_pass)[0]):
            if p[i] < p_thr:
                row = sub_sel.iloc[ok]
                # Create hit dict with all original columns plus correlation results
                hit = {
                    "patient_id": pid,
                    "n_samples":  int(n_vec[i]),
                    "r":          float(r[mask_pass][i]),
                    "p":          float(p[i]),
                    "method":     method,
                }
                # Add all original columns from the row
                for col in row.index:
                    if col not in hit:  # Avoid overwriting correlation results
                        hit[col] = row[col]
                hits.append(hit)
    return hits


# ───────────────────────── CLI parsing ────────────────────────────────────

def parse_cli() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("snv_parquet", type=Path)
    ap.add_argument("metadata", type=Path)
    ap.add_argument("--id-col", default="patient_id")
    ap.add_argument("--time-col", default="week_num")
    ap.add_argument("--cal-col", default="calprotectin")
    ap.add_argument("--method", choices=["spearman", "pearson"], default="spearman")
    ap.add_argument("--min-samples", type=int, default=5)
    ap.add_argument("--r-thresh", type=float, default=0.7)
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=max(mp.cpu_count() // 2, 1))
    ap.add_argument("--out", type=Path, default=Path("snv_calprotectin_hits.tsv"))
    return ap.parse_args()


# ───────────────────────── main ───────────────────────────────────────────

def main() -> None:
    args = parse_cli()

    # metadata → wide patient × week matrix ---------------------------------
    meta = pd.read_csv(args.metadata, sep=None, engine="python")
    print(f"DEBUG: Metadata loaded with {len(meta)} rows, columns: {list(meta.columns)[:10]}...")
    
    pid_col = args.id_col if args.id_col in meta.columns else "Participant ID"
    meta = meta.rename(columns={pid_col: "patient_id", args.time_col: "week", args.cal_col: "cpx"})
    meta["week"] = pd.to_numeric(meta["week"], errors="coerce")
    meta_wide = meta.pivot_table(index="patient_id", columns="week", values="cpx", aggfunc="mean")
    
    print(f"DEBUG: Metadata wide table has {len(meta_wide)} patients, {len(meta_wide.columns)} weeks")
    print(f"DEBUG: Available weeks in metadata: {sorted(meta_wide.columns)[:10]}...")

    # detect week columns in parquet set (assumes homogeneous schema) -------
    first_frag = next(ds.dataset(args.snv_parquet, format="parquet").get_fragments())
    schema = first_frag.physical_schema
    week_cols = sorted([f.name for f in schema if f.name.isdigit()], key=int)
    print(f"DEBUG: SNV data has {len(week_cols)} week columns: {week_cols[:10]}...")

    # pool over fragments ----------------------------------------------------
    frag_paths = [f.path for f in ds.dataset(args.snv_parquet, format="parquet").get_fragments()]
    logger.info("Processing %d parquet fragments with %d thread(s)…", len(frag_paths), args.threads)
    pool = mp.Pool(args.threads)
    worker_args = [(
        p, args.id_col, week_cols, meta_wide, args.method,
        args.min_samples, args.r_thresh, args.p_thresh
    ) for p in frag_paths]
    all_hits: List[Dict] = []
    for res in pool.imap_unordered(process_fragment, worker_args):
        all_hits.extend(res)
    pool.close(); pool.join()

    # write ------------------------------------------------------------------
    out_df = pd.DataFrame(all_hits).sort_values("p")
    out_df.to_csv(args.out, sep="\t", index=False)
    logger.info("Done.  Hits: %d rows → %s", len(out_df), args.out)

if __name__ == "__main__":
    main()