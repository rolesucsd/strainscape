#!/usr/bin/env python3
"""
Prepare mutation trajectories for visualisation (fast version).

Inputs:
  – mutations.tsv  (per-sample counts)
  – trends.tsv     (per-site model fits)
  – metadata.csv   (Sample ↔ week_num)
Outputs
  mutation_trajectories.tsv   (wide table)
  mutation_trends.tsv         (verbatim copy of trends)
"""

from pathlib import Path
import numpy as np, pandas as pd, argparse, logging, sys, os
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

# ─────────── main routine ────────────
def prepare_trajectories_fast(muts_f, trends_f, meta_f, out_dir):
    logger.info("Reading tables …")
    meta   = (pd.read_csv(meta_f, usecols=["External.ID", "week_num"])
                .rename(columns={"External.ID": "Sample"}))

    trends = pd.read_csv(trends_f, sep="\t")

    muts   = pd.read_csv(muts_f,   sep="\t")

    # ── restrict to sites present in trends ──
    muts = muts.merge(
        trends[["scaffold", "position"]].drop_duplicates(),
        on=["scaffold", "position"],
        how="inner",
    )

    # ── merge metadata ──
    muts = (muts.merge(meta,  on="Sample", how="left")
                  .dropna(subset=["week_num"]))
    muts["week_num"] = muts["week_num"].astype(int)

    # ── compute alt-allele freq vectorised ──
    base_cols = ["A", "C", "G", "T"]
    base_to_idx = dict(zip(base_cols, range(len(base_cols))))

    counts = muts[base_cols].to_numpy(dtype=float)
    idx    = muts["new_base"].map(base_to_idx).to_numpy(int, na_value=-1)

    freq = np.full(len(muts), np.nan, dtype=float)
    valid = idx >= 0
    freq[valid] = counts[np.arange(len(muts))[valid], idx[valid]] / muts.loc[valid, "position_coverage"]
    muts["Value"] = freq

    # ── reshape to wide trajectory table ──
    key_cols = ["scaffold", "position", "ref_base", "new_base"]
    traj = (muts
            .loc[:, key_cols + ["week_num", "Value"]]
            .dropna(subset=["Value"])
            .pivot_table(index=key_cols,
                         columns="week_num",
                         values="Value",
                         aggfunc="first")
            .reset_index())

    # add mutation type & model fits
    traj["Type"] = traj["ref_base"] + ">" + traj["new_base"]
    traj = traj.merge(
        trends.rename(columns={"slope": "OLS_slope", "p_value": "OLS_pvalue",
                               "r_squared": "OLS_fit"}),
        on=["scaffold", "position"],
        how="left",
    )

    # final tidy / column order
    numeric_cols = sorted(c for c in traj.columns if isinstance(c, (int, np.integer)))
    ordered = (["scaffold", "position", "ref_base", "new_base", "Type",
                "OLS_slope", "OLS_pvalue", "OLS_fit"] + numeric_cols)
    traj = traj[[c for c in ordered if c in traj]]

    # ── write outputs ──
    out_dir.mkdir(parents=True, exist_ok=True)
    traj.to_csv(out_dir / "mutation_trajectories.tsv", sep="\t", index=False)
    trends.to_csv(out_dir / "mutation_trends.tsv",       sep="\t", index=False)

    logger.info(f"Trajectories: {len(traj):,} rows → {out_dir/'mutation_trajectories.tsv'}")


# ───────────── CLI ─────────────
if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--mutation_file", required=True)
    p.add_argument("--trends_file",    required=True)
    p.add_argument("--metadata_file",  required=True)
    p.add_argument("--output_dir",     required=True)
    p.add_argument("--log_file")
    args = p.parse_args()

    if args.log_file:
        setup_logging(args.log_file)

    prepare_trajectories_fast(Path(args.mutation_file),
                              Path(args.trends_file),
                              Path(args.metadata_file),
                              Path(args.output_dir))
