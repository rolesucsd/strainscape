#!/usr/bin/env python3
"""hotspot_finder.py

Identify mutation hotspots in SNV parquet tables.

---------------------------------------------------------------------
Assumptions about the parquet schema
---------------------------------------------------------------------
Your SNV dataset (single parquet directory) must contain **at least**
these columns:
    patient_id      – optional, ignored for hotspot calling
    gene            – gene identifier (string)
    gene_length     – length of that gene in bp (int)
    strand          – '+' or '-' (string)
    position        – position of SNV *inside* the gene (0‑based, int)
    dbxrefss          – dbxrefss / locus_tag / whatever identifier you need
If your column names differ, map them with CLI flags (see --help).

---------------------------------------------------------------------
Outputs (always written next to the parquet directory):
---------------------------------------------------------------------
1.  <out_prefix>_hotspots.txt
    Columns: gene, gene_length, strand, hotspot_start, hotspot_end, count, dbxrefs

2.  <out_prefix>_window_counts.txt
    All 10‑bp windows for every gene, with mutation counts – tidy format:
        gene, window_start, window_end, count, gene_length, dbxrefs
    Ready to plot as a bar‑plot / heatmap in R.

---------------------------------------------------------------------
Hotspot definition
---------------------------------------------------------------------
Within each gene, we slide a **fixed window (default 10 bp)** across the
strand‑normalised coordinate vector v = [0 … gene_length‑1].

For a given gene we call a window a *hotspot* when its mutation count is
>= `quantile_cutoff` (default 0.99) of all window counts in that gene.
Adjust with --cutoff or --absolute if you prefer a fixed threshold.

---------------------------------------------------------------------
Usage examples
---------------------------------------------------------------------
$ python hotspot_finder.py /path/to/all_snvs \
        --gene-col gene --pos-col position --length-col gene_length \
        --strand-col strand --dbxrefs-col dbxrefs \
        --window 10 --cutoff 0.99

$ python hotspot_finder.py all_snvs --absolute 5   # any window >=5 muts

---------------------------------------------------------------------
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────────────
# helper
# ────────────────────────────────────────────────────────────────────

def strand_normalise(pos: np.ndarray, gene_length: int, strand: str) -> np.ndarray:
    """Return positions always counting 0→gene_length‑1 from 5' to 3'."""
    if strand == "-":
        return gene_length - 1 - pos
    return pos


def window_counts(positions: np.ndarray, gene_length: int, window: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return (bin_edges, counts) for a fixed‑width histogram."""
    bins = np.arange(0, gene_length + window, window)
    counts, _ = np.histogram(positions, bins=bins)
    return bins, counts


# ────────────────────────────────────────────────────────────────────
# main logic
# ────────────────────────────────────────────────────────────────────

def process_gene(df: pd.DataFrame, window: int, cutoff_q: float | None, absolute: int | None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    gene = df.iloc[0]["gene"]
    gene_length = int(df.iloc[0]["gene_length"])
    strand = df.iloc[0]["strand"]
    dbxrefs = df.iloc[0]["dbxrefs"]

    positions = df["position"].to_numpy(dtype=int)
    norm_pos = strand_normalise(positions, gene_length, strand)

    bins, counts = window_counts(norm_pos, gene_length, window)

    # tidy long‑format for plotting
    window_df = pd.DataFrame({
        "gene": gene,
        "window_start": bins[:-1],
        "window_end": bins[1:],
        "count": counts,
        "gene_length": gene_length,
        "dbxrefs": dbxrefs,
    })

    # hotspot cut‑off ------------------------------
    if absolute is not None:
        mask = counts >= absolute
    else:
        threshold = np.quantile(counts, cutoff_q) if counts.size else 0
        mask = counts >= threshold

    hotspot_df = window_df[mask].copy()
    hotspot_df.rename(columns={"window_start": "hotspot_start", "window_end": "hotspot_end"}, inplace=True)
    return hotspot_df, window_df


# ────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Detect mutation hotspots in an SNV parquet dataset")
    p.add_argument("parquet", type=Path, help="Path to SNV parquet directory")
    p.add_argument("--window", type=int, default=10, help="Sliding window size (bp)")
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--cutoff", type=float, default=0.99, help="Quantile cutoff (e.g. 0.99) within gene")
    grp.add_argument("--absolute", type=int, default=5, help="Absolute count threshold (e.g. 5 mutations)")
    p.add_argument("--gene-col", default="gene")
    p.add_argument("--pos-col", default="position")
    p.add_argument("--start", default="start")
    p.add_argument("--stop", default="stop")
    p.add_argument("--strand-col", default="strand")
    p.add_argument("--dbxrefs-col", default="dbxrefs")
    p.add_argument("--out-prefix", default="hotspots", help="Prefix for output files")
    return p.parse_args()


def main():
    args = parse_args()

    # read parquet -------------------
    logger.info("Loading parquet from %s", args.parquet)
    dset = ds.dataset(args.parquet, format="parquet")
    tbl = dset.to_table()
    df = tbl.to_pandas()

    df = df[df['coding_status'] == "Coding"]

    col_map: Dict[str, str] = {
        "gene": args.gene_col,
        "position": args.pos_col,
        "start": args.start,
        "stop": args.stop,
        "strand": args.strand_col,
        "dbxrefs": args.dbxrefs_col,
    }
    df = df.rename(columns={v: k for k, v in col_map.items()})
    if "start" in df.columns and "stop" in df.columns:
        df["gene_length"] = df["stop"].astype(float) - df["start"].astype(float)
        logger.info("Computed gene_length from start/stop columns")

    required = ["gene", "position", "gene_length", "strand", "dbxrefs"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {', '.join(missing)}")

    logger.info("Processing %d SNVs across %d genes", len(df), df["gene"].nunique())

    hotspot_list: List[pd.DataFrame] = []
    window_list: List[pd.DataFrame] = []

    for (gene, glen), gdf in df.groupby(["gene", "gene_length"], sort=False):
        h, w = process_gene(gdf, window=args.window, cutoff_q=args.cutoff, absolute=args.absolute)
        hotspot_list.append(h)
        window_list.append(w)

    hotspot_df = pd.concat(hotspot_list, ignore_index=True) if hotspot_list else pd.DataFrame()
    window_df = pd.concat(window_list, ignore_index=True) if window_list else pd.DataFrame()

    # write outputs ------------------
    prefix = Path(args.out_prefix)
    hotspot_path = Path(f"{prefix}_hotspots.txt")
    window_path  = Path(f"{prefix}_window_counts.txt")
    hotspot_df.to_csv(hotspot_path, sep="\t", index=False)
    window_df.to_csv(window_path, sep="\t", index=False)

    logger.info("Hotspots → %s (%d rows)", hotspot_path, len(hotspot_df))
    logger.info("Window counts → %s (%d rows)", window_path, len(window_df))


if __name__ == "__main__":
    main()
