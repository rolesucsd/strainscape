#!/usr/bin/env python3
"""
Classify evolutionary events into SLAM / CGR / SCT (primary classes)
from ALL mutations (no high-AF pre-filter), then overlay existing
high-ΔAF seeds from analyzed_mutations.py for reporting.

Inputs (TSV; column names shown are expected; adapt if your schema differs):
  snv.tsv       : [bin_id, scaffold, position, time, af, depth]
  linkage.tsv   : [bin_id, scaffold, pos1, pos2, r2, time]
  sv.tsv        : [bin_id, scaffold, sv_start, sv_end, sv_type, support_split]
  analyzed_mutations.tsv (optional): [bin_id, scaffold, position, gene_label, delta_af, ...]
Outputs:
  blocks.tsv            : time series (time, block_id, af)
  events.tsv            : per-block metadata + class
  seeds_with_class.tsv  : your seeds mapped to block & class
"""

from __future__ import annotations
import pandas as pd, numpy as np, argparse, json, sys
from pathlib import Path

# -------- thresholds (pre-register these) ----------
R2_HARD       = 0.90     # linkage threshold
MAX_LINK_DIST = 10_000   # bp
MIN_DEPTH     = 20
MIN_TPTS      = 3        # min timepoints per block to keep trajectory
CGR_MIN_BP    = 3_000    # size-only CGR criterion (or use n_snvs)
CGR_MIN_SNVS  = 5
SV_SPLIT_MIN  = 5
SCT_PC1_MIN_VAR = 0.40   # genome-wide coherence
SCT_LOADING_MIN = 0.20   # absolute PC1 loading to label a block as part of SCT

# ---------------- I/O ----------------
def load_table(path: str) -> pd.DataFrame:
    if not path: return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str)
    # numeric coercions (best-effort)
    for c in ("position","pos1","pos2","sv_start","sv_end","time","af","depth","r2","support_split"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

# ------------- block building -------------
def build_ld_blocks_all_times(link_all: pd.DataFrame) -> dict[tuple,str]:
    """
    Build LD blocks per (bin_id, scaffold) using *union* of strong/near edges across times.
    Returns dict key=(bin,scaf) -> list of blocks, each as sorted list of positions.
    """
    blocks_by_key: dict[tuple,list] = {}
    if link_all.empty: return blocks_by_key
    # keep strong/near edges
    L = link_all[(link_all["r2"] >= R2_HARD) &
                 (link_all["pos1"] - link_all["pos2"]).abs().le(MAX_LINK_DIST)].copy()
    for (bin_id, scaf), g in L.groupby(["bin_id","scaffold"], dropna=False):
        # adjacency lists
        adj: dict[int,set] = {}
        for r in g.itertuples():
            p1, p2 = int(r.pos1), int(r.pos2)
            adj.setdefault(p1,set()).add(p2)
            adj.setdefault(p2,set()).add(p1)
        nodes = set(adj.keys())
        comps = []
        while nodes:
            start = nodes.pop()
            stack = [start]; comp = {start}
            while stack:
                u = stack.pop()
                for v in adj.get(u,()):
                    if v in nodes:
                        nodes.remove(v); comp.add(v); stack.append(v)
            comps.append(sorted(comp))
        blocks_by_key[(bin_id, scaf)] = comps
    return blocks_by_key

def summarize_block_timeseries(snv: pd.DataFrame, bin_id: str, scaf: str, positions: list[int]):
    df = snv[(snv["bin_id"]==bin_id) & (snv["scaffold"]==scaf) &
             (snv["position"].isin(positions))].copy()
    if df.empty: return None, None, None
    df = df[df["depth"] >= MIN_DEPTH]
    if df.empty: return None, None, None
    # depth-weighted AF per time
    ts = (df.groupby("time")
            .apply(lambda g: np.average(g["af"], weights=g["depth"]))
            .rename("af").reset_index())
    if ts["time"].nunique() < MIN_TPTS: return None, None, None
    x, y = ts["time"].values, ts["af"].values
    slope = np.polyfit(x, y, 1)[0] if len(ts)>=2 else 0.0
    delta = float(np.nanmax(y) - np.nanmin(y))
    return ts.set_index("time")["af"], slope, delta

def overlap(a0,a1,b0,b1): return (a0<=b1) and (b0<=a1)

def sv_tags_for_span(sv: pd.DataFrame, bin_id: str, scaf: str, s0: int, s1: int) -> list[str]:
    tags = []
    if sv.empty: return tags
    H = sv[(sv["bin_id"]==bin_id) & (sv["scaffold"]==scaf) &
           sv.apply(lambda r: overlap(s0,s1, r["sv_start"], r["sv_end"]), axis=1)]
    if not H.empty:
        tags.append("SV")
        for t in H["sv_type"].dropna().unique():
            tags.append(f"SV_{t}")
        if (H["support_split"] >= SV_SPLIT_MIN).any(): tags.append("SV_split_support")
    return sorted(set(tags))

def classify_primary(n_snvs, span_bp, has_sv):
    if has_sv or n_snvs >= CGR_MIN_SNVS or span_bp >= CGR_MIN_BP:
        return "CGR"
    elif n_snvs <= 2 and span_bp <= 2000:
        return "SLAM"
    else:
        return "local_multilocus"

# ------------- SCT detection (per bin) -------------
def pca_first_component(X: pd.DataFrame):
    # standardize columns
    Z = (X - X.mean())/X.std(ddof=0)
    Z = Z.replace([np.inf,-np.inf], np.nan).fillna(0.0)
    U,S,Vt = np.linalg.svd(Z.values, full_matrices=False)
    var_expl = (S[0]**2) / (S**2).sum() if S.sum()>0 else 0.0
    pc1_time = pd.Series(U[:,0]*S[0], index=Z.index)
    loadings = pd.Series(Vt[0,:], index=Z.columns)
    return float(var_expl), pc1_time, loadings

# ------------- main driver -------------
def run(snv_path, link_path, sv_path, seeds_path, outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    snv = load_table(snv_path)
    link = load_table(link_path)
    sv   = load_table(sv_path) if sv_path else pd.DataFrame()
    seeds = load_table(seeds_path) if seeds_path else pd.DataFrame()

    required = [["bin_id","scaffold","position","time","af","depth"],
                ["bin_id","scaffold","pos1","pos2","r2"]]
    for req, df, nm in zip(required, [snv, link], ["snv","linkage"]):
        miss = [c for c in req if c not in df.columns]
        if miss:
            sys.exit(f"[error] {nm} is missing columns: {miss}")

    # Build LD blocks (union across times)
    blocks_union = build_ld_blocks_all_times(link)

    # Assemble per-block metadata & time series
    block_rows = []
    block_ts_cols = {}  # block_id -> pd.Series(time->af)
    for (bin_id, scaf), comps in blocks_union.items():
        # fallback: treat singletons (no edges) that appear only in SNV but not linkage
        if not comps:
            # add singleton positions present in SNV
            positions = snv[(snv["bin_id"]==bin_id)&(snv["scaffold"]==scaf)]["position"].dropna().astype(int).unique()
            comps = [[int(p)] for p in positions]
        for i, pos_list in enumerate(comps, start=1):
            pos_list = sorted(set(int(p) for p in pos_list))
            s0, s1 = (min(pos_list), max(pos_list))
            series, slope, delta = summarize_block_timeseries(snv, bin_id, scaf, pos_list)
            if series is None: 
                continue
            block_id = f"{bin_id}|{scaf}|{s0}-{s1}"
            block_ts_cols[block_id] = series
            tags = sv_tags_for_span(sv, bin_id, scaf, s0, s1)
            primary = classify_primary(n_snvs=len(pos_list), span_bp=s1-s0, has_sv=("SV" in tags))
            block_rows.append({
                "block_id": block_id, "bin_id": bin_id, "scaffold": scaf,
                "span_start": s0, "span_end": s1, "span_bp": s1-s0,
                "n_snvs": len(pos_list), "primary_class_preSCT": primary,
                "sv_tags": ",".join(tags),
                "slope": slope, "delta_af": delta
            })

    if not block_rows:
        sys.exit("[error] no blocks with sufficient timepoints after depth filtering")

    events = pd.DataFrame(block_rows)

    # Build per-bin matrices and detect SCT; tag blocks with SCT if they contribute to PC1
    events["primary_class"] = events["primary_class_preSCT"]
    blocks_long = []
    for bin_id, ebin in events.groupby("bin_id"):
        # time × blocks matrix
        cols = {bid: block_ts_cols[bid] for bid in ebin["block_id"]}
        mat = pd.DataFrame(cols).sort_index()
        if mat.shape[1] < 3 or mat.shape[0] < MIN_TPTS:
            # not enough blocks or timepoints to evaluate SCT
            pass
        else:
            var_expl, pc1_time, load = pca_first_component(mat)
            if var_expl >= SCT_PC1_MIN_VAR:
                # select blocks with strong loadings (same sign majority automatically captured)
                sel = load[load.abs() >= SCT_LOADING_MIN].index
                events.loc[events["block_id"].isin(sel) & (events["bin_id"]==bin_id), "primary_class"] = "SCT"
        # collect long trajectories
        long = mat.stack().rename("af").reset_index().rename(columns={"level_0":"time","level_1":"block_id"})
        long["bin_id"] = bin_id
        blocks_long.append(long)

    blocks_long = pd.concat(blocks_long, ignore_index=True)

    # Overlay your ΔAF seeds (if provided): map each seed position into a block span
    if not seeds.empty:
        seeds["position"] = pd.to_numeric(seeds["position"], errors="coerce")
        merged = seeds.merge(
            events[["bin_id","scaffold","block_id","span_start","span_end","primary_class"]],
            on=["bin_id","scaffold"], how="left"
        )
        in_block = (merged["position"] >= merged["span_start"]) & (merged["position"] <= merged["span_end"])
        seeds_with_class = merged[in_block].drop(columns=["span_start","span_end"])
        seeds_with_class.to_csv(outdir/"seeds_with_class.tsv", sep="\t", index=False)

    # Write outputs
    events.to_csv(outdir/"events.tsv", sep="\t", index=False)
    blocks_long.to_csv(outdir/"blocks.tsv", sep="\t", index=False)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Classify evolutionary events into SLAM / CGR / SCT.")
    ap.add_argument("--snv", required=True)
    ap.add_argument("--linkage", required=True)
    ap.add_argument("--sv", required=False)
    ap.add_argument("--seeds", required=False, help="Output from analyzed_mutations.py (ΔAF seeds)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    run(args.snv, args.linkage, args.sv, args.seeds, args.outdir)
