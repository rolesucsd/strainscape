#!/usr/bin/env python3
# calculate_combined_stats_v9.py  ───────────────────────────────────────────
#
# • Streams a huge mutation table in chunks
# • *Keeps only sweeping SNVs*   (min_freq < 0·2  &  max_freq > 0·8)
# • Writes one big Parquet with full annotation for those sweeps
# • Detects sweep variants based on frequency changes over time
# ---------------------------------------------------------------------------

import argparse, gzip, re, uuid, sys
from pathlib import Path
from typing   import Dict, List, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds
from tqdm import tqdm

# ── CLI ────────────────────────────────────────────────────────────────────
def get_args():
    p = argparse.ArgumentParser("iHMP sweeps detection streamer")
    p.add_argument("--mut-file",      required=True,
                   help="combined_trend_mapped_type.txt[.gz]")
    p.add_argument("--scaffold-file", required=True,
                   help="combined_processed_scaffolds.txt")
    p.add_argument("--meta-file",     required=True,
                   help="meta_map.csv  (only to map sample→patient)")
    p.add_argument("--out-dir",       default="mutation_to_snv_sweep_parquet")
    p.add_argument("--chunk-size",    type=int, default=1_000_000)
    p.add_argument("--low-mem",       action="store_true",
                   help="pass low_memory=False to pandas read_csv")
    return p.parse_args()

# ── helpers ────────────────────────────────────────────────────────────────
def detect_week_cols(cols: List[str]) -> List[int]:
    return sorted(int(c) for c in cols if re.fullmatch(r"\d+", c))


def split_dbxrefs(s: str):
    return [t.strip().lower() for t in str(s).split(",") if t.strip()]

def extract_go(tok):   return [t for t in tok if t.startswith("go:")]
def extract_ec(tok): return [t for t in tok if t.startswith("ec:") or re.fullmatch(r"k\d{5}", t)]
def extract_kegg(tok): return [t for t in tok if t.startswith("kegg:") or re.fullmatch(r"k\d{5}", t)]

def load_bin_map(scaffold_file: str, meta_file: str) -> Dict[Tuple[str,str], str]:
    print("[DEBUG] Reading meta file:", meta_file)
    meta = pd.read_csv(meta_file, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(" ","_")

    # Prefer 'run' if present (per updated metadata), else fall back to External ID variants
    candidate_sample_cols = ["run", "external.id", "external_id"]
    sample_col = None
    for c in candidate_sample_cols:
        if c in meta.columns and meta[c].notna().any():
            sample_col = c
            break
    if sample_col is None:
        raise ValueError("Could not find a sample identifier column in metadata (tried: Run, External.ID).")

    # Harmonise columns and drop NAs to avoid merge issues
    meta = (meta.rename(columns={sample_col: "sample"})
                [["sample", "participant_id"]]
                .rename(columns={"participant_id": "patient_id"})
                .dropna(subset=["sample", "patient_id"]))

    # Normalise whitespace/case for robust merges
    meta["sample"] = meta["sample"].astype(str).str.strip()
    meta["patient_id"] = meta["patient_id"].astype(str).str.strip()
    print(f"[DEBUG] Meta rows: {len(meta)} | unique samples: {meta['sample'].nunique()} | unique patients: {meta['patient_id'].nunique()}")

    print("[DEBUG] Reading scaffold file:", scaffold_file)
    scaf = pd.read_csv(scaffold_file, sep="\t",
                       usecols=["scaffold","Sample","bin"], low_memory=False)
    scaf.columns = scaf.columns.str.lower()

    # Ensure sample column normalisation matches metadata
    if "sample" in scaf.columns:
        scaf["sample"] = scaf["sample"].astype(str).str.strip()
    print(f"[DEBUG] Scaffold rows: {len(scaf)} | unique scaffolds: {scaf['scaffold'].nunique()} | unique samples: {scaf['sample'].nunique()} | unique bins: {scaf['bin'].nunique()}")

    print("[DEBUG] Merging scaffold↔meta on sample")
    mp_full = scaf.merge(meta, on="sample", how="left")
    missing_pat = mp_full["patient_id"].isna().sum()
    print(f"[DEBUG] Merge rows: {len(mp_full)} | rows missing patient_id: {missing_pat}")
    mp = mp_full.dropna(subset=["patient_id"])
    print(f"[DEBUG] Mapping entries (patient, scaffold)->bin: {len(mp)} (unique keys: {mp[['patient_id','scaffold']].drop_duplicates().shape[0]})")
    return {(r.patient_id, r.scaffold): r.bin for r in mp.itertuples(index=False)}

def append_parquet(dir_path: Path, df: pd.DataFrame):
    if df.empty: return
    dir_path.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   dir_path/f"part-{uuid.uuid4().hex}.parquet",
                   compression="snappy")

# ── streaming ──────────────────────────────────────────────────────────────
def main():
    a = get_args()
    out_dir = Path(a.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # ── prep ----------------------------------------------------------------
    preview   = pd.read_csv(a.mut_file, sep="\t", nrows=2000, low_memory=False)
    week_cols = detect_week_cols(preview.columns)
    if not week_cols:
        sys.exit("❌ No numeric week columns (weeks) detected.")
    print(f"Detected {len(week_cols)} week columns → {week_cols[:10]} …")

    bin_map = load_bin_map(a.scaffold_file, a.meta_file)

    opener      = gzip.open if a.mut_file.endswith(".gz") else open
    total_rows  = sum(1 for _ in opener(a.mut_file, "rt")) - 1
    pbar        = tqdm(total=total_rows, unit="rows")


    # ── iterate CSV chunks ---------------------------------------------------
    for chunk in pd.read_csv(a.mut_file, sep="\t",
                             chunksize=a.chunk_size,
                             low_memory=not a.low_mem):
        chunk.columns = (chunk.columns.str.strip()
                                       .str.lower()
                                       .str.replace(" ","_"))

        # minimal numeric coercion ------------------------------------------
        num_cols = ["mutation_type","ref_base","new_base",
                    *map(str,week_cols)]
        for c in map(str, week_cols):
            if c in chunk.columns:
                chunk[c] = pd.to_numeric(chunk[c], errors="coerce")

        # annotate bin -------------------------------------------------------
        if "bin" not in chunk.columns:
            chunk["bin"] = [bin_map.get((pid, chrom))
                            for pid, chrom in zip(chunk.patient_id, chunk.chromosome)]
        # Chunk-level debug on bin assignment
        na_bins = chunk["bin"].isna().sum() if "bin" in chunk.columns else len(chunk)
        uniq_bins = chunk["bin"].nunique(dropna=True) if "bin" in chunk.columns else 0
        print(f"[DEBUG] Chunk rows: {len(chunk)} | NA bins: {na_bins} | unique bins in chunk: {uniq_bins}")

        # functional tokens --------------------------------------------------
        if "dbxrefs" in chunk.columns and "go_terms" not in chunk.columns:
            toks = chunk.dbxrefs.apply(split_dbxrefs)
            chunk["go_terms"]   = toks.apply(extract_go)
            chunk["kegg_terms"] = toks.apply(extract_kegg)
            chunk["ec_terms"] = toks.apply(extract_ec)

        # per-variant min/max freq ------------------------------------------
        # 1) Make sure week columns exist as strings and are numeric
        week_str = [str(c) for c in week_cols if str(c) in chunk.columns]
        for c in week_str:
            chunk[c] = pd.to_numeric(chunk[c], errors="coerce")

        # 2) Auto-detect percent scale and clamp to [0,1]
        vals = chunk[week_str].to_numpy(dtype=float)  # shape = (n_rows, n_weeks)

        # Heuristic: if ≥10% of entries > 1 and <1% > 100, treat as percent & scale
        with np.errstate(invalid="ignore"):
            frac_gt1   = np.nanmean(vals > 1.0)
            frac_gt100 = np.nanmean(vals > 100.0)
        if (frac_gt1 > 0.10) and (frac_gt100 < 0.01):
            vals = vals / 100.0

        # Clamp to [0,1]
        vals = np.clip(vals, 0.0, 1.0)

        # 3) Per-variant stats across weeks
        n_obs         = np.sum(~np.isnan(vals), axis=1)
        min_f         = np.nanmin(vals, axis=1)
        max_f         = np.nanmax(vals, axis=1)
        freq_range    = max_f - min_f

        # 4) Put back into the chunk (after clamping)
        chunk["n_obs"]      = n_obs
        chunk["min_f"]      = min_f
        chunk["max_f"]      = max_f
        chunk["freq_range"] = freq_range

        # 5) Robust sweep call (direction-agnostic)
        LO, HI        = 0.20, 0.80
        MIN_POINTS    = 3         # require at least 3 observed time points (was 5, but 27.5% of variants only have 3 obs)
        MIN_RANGE     = 0.60      # HI-LO is 0.60, so anything that meets extremes will pass

        sweep_mask = (chunk["n_obs"] >= MIN_POINTS) & \
                     (chunk["min_f"] <= LO) & (chunk["max_f"] >= HI) & \
                     (chunk["freq_range"] >= MIN_RANGE)

        chunk["is_sweep"] = sweep_mask

        append_parquet(out_dir / "all_snvs",   chunk)
        append_parquet(out_dir / "sweep_variants", chunk.loc[sweep_mask])

        pbar.update(len(chunk))
    pbar.close()

    # ────────────────────────────────────────────────────────────────────────
    print("🏁 Finished – wrote:")
    print("   • all_snvs/part-*.parquet        (every SNV, flag = is_sweep)")
    print("   • sweep_variants/part-*.parquet  (subset where is_sweep == TRUE)")

if __name__ == "__main__":
    main()
