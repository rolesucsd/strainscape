#!/usr/bin/env python3
# calculate_combined_stats_v9.py  ───────────────────────────────────────────
#
# • Streams a huge mutation table in chunks
# • *Keeps only sweeping SNVs*   (min_freq < 0·2  &  max_freq > 0·8)
# • Writes one big Parquet with full annotation for those sweeps
# • Writes per-gene and per-bin dN/dS tables (all coding SNVs)
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
    p = argparse.ArgumentParser("iHMP sweeps + dN/dS streamer")
    p.add_argument("--mut-file",      required=True,
                   help="combined_trend_mapped_type.txt[.gz]")
    p.add_argument("--scaffold-file", required=True,
                   help="combined_processed_scaffolds.txt")
    p.add_argument("--meta-file",     required=True,
                   help="hmp2_metadata_2018-08-20.csv  (only to map sample→patient)")
    p.add_argument("--out-dir",       default="mutation_to_snv_sweep_parquet")
    p.add_argument("--chunk-size",    type=int, default=1_000_000)
    p.add_argument("--low-mem",       action="store_true",
                   help="pass low_memory=False to pandas read_csv")
    return p.parse_args()

# ── helpers ────────────────────────────────────────────────────────────────
def detect_week_cols(cols: List[str]) -> List[int]:
    return sorted(int(c) for c in cols if re.fullmatch(r"\d+", c))

def classify_dn_ds(mut_type: pd.Series):
    m = mut_type.str.lower().fillna("")
    return m.isin(["missense", "nonsense"]), m.eq("silent")

def split_dbxrefs(s: str):
    return [t.strip().lower() for t in str(s).split(",") if t.strip()]

def extract_go(tok):   return [t for t in tok if t.startswith("go:")]
def extract_ec(tok): return [t for t in tok if t.startswith("ec:") or re.fullmatch(r"k\d{5}", t)]
def extract_kegg(tok): return [t for t in tok if t.startswith("kegg:") or re.fullmatch(r"k\d{5}", t)]

def load_bin_map(scaffold_file: str, meta_file: str) -> Dict[Tuple[str,str], str]:
    meta = pd.read_csv(meta_file, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(" ","_")
    sample_col  = "external.id" if "external.id" in meta.columns else "external_id"
    meta = meta.rename(columns={sample_col:"sample"})[["sample","participant_id"]]\
               .rename(columns={"participant_id":"patient_id"}).dropna()

    scaf = pd.read_csv(scaffold_file, sep="\t",
                       usecols=["scaffold","Sample","bin"], low_memory=False)
    scaf.columns = scaf.columns.str.lower()

    mp = scaf.merge(meta, on="sample", how="left").dropna(subset=["patient_id"])
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

    runs = {"gene_dnds":None, "bin_dnds":None}       # incremental dN/dS

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

        # functional tokens --------------------------------------------------
        if "dbxrefs" in chunk.columns and "go_terms" not in chunk.columns:
            toks = chunk.dbxrefs.apply(split_dbxrefs)
            chunk["go_terms"]   = toks.apply(extract_go)
            chunk["kegg_terms"] = toks.apply(extract_kegg)
            chunk["ec_terms"] = toks.apply(extract_ec)

        # per-variant min/max freq ------------------------------------------
        id_cols = [c for c in chunk.columns if c not in map(str,week_cols)]
        long_df = (chunk.melt(id_vars=id_cols, value_vars=list(map(str, week_cols)),
                              var_name="week_num", value_name="freq")
                         .dropna(subset=["freq"]))
        var_mm = (long_df.groupby(["patient_id","chromosome","position"], as_index=False)
                           .agg(min_f=("freq","min"),
                                max_f=("freq","max")))

        # join min/max back to main table  ----------------------------------
        chunk = chunk.merge(var_mm, on=["patient_id","chromosome","position"])

        freq_cols = list(map(str, week_cols))           # ensure week columns are strings
        # apply row-wise; fine inside a 1 M-row chunk
        sweep_mask            = (chunk.min_f < .20) & (chunk.max_f > .80)
        chunk["is_sweep"]     = sweep_mask       # boolean column

        append_parquet(out_dir / "all_snvs",   chunk)
        append_parquet(out_dir / "sweep_variants", chunk.loc[sweep_mask])

        # dN/dS – still use ALL coding SNVs ---------------------------------
        is_dn, is_ds = classify_dn_ds(chunk.mutation_type)
        coding_mask  = is_dn | is_ds
        coding       = chunk[coding_mask]

        if "gene" in coding.columns:
            gene_part = (coding.assign(is_dn=is_dn[coding_mask],
                                       is_ds=is_ds[coding_mask])
                               .groupby(["patient_id","gene"], as_index=False)
                               .agg(dn=("is_dn","sum"), ds=("is_ds","sum")))
            runs["gene_dnds"] = (pd.concat([runs["gene_dnds"], gene_part])
                                   if runs["gene_dnds"] is not None else gene_part)

        if "bin" in coding.columns:
            bin_part  = (coding.assign(is_dn=is_dn[coding_mask],
                                       is_ds=is_ds[coding_mask])
                               .dropna(subset=["bin"])
                               .groupby(["patient_id","bin"], as_index=False)
                               .agg(dn=("is_dn","sum"), ds=("is_ds","sum")))
            runs["bin_dnds"] = (pd.concat([runs["bin_dnds"], bin_part])
                                   if runs["bin_dnds"] is not None else bin_part)

        pbar.update(len(chunk))
    pbar.close()

    # ── finalise dN/dS tables ----------------------------------------------
    for name in ["gene_dnds","bin_dnds"]:
        df = runs[name]
        if df is None: continue
        df = df.groupby(df.columns.tolist()[:-2], as_index=False).sum()  # collapse duplicates
        df["dnds"] = df.apply(lambda r: r.dn / r.ds if r.ds else np.nan, axis=1)
        pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                       out_dir/f"{name}.parquet", compression="snappy")
        print(f"✔ {name}.parquet ({df.shape[0]:,} rows)")

    # ────────────────────────────────────────────────────────────────────────
    print("🏁 Finished – wrote:")
    print("   • all_snvs/part-*.parquet        (every SNV, flag = is_sweep)")
    print("   • sweep_variants/part-*.parquet  (subset where is_sweep == TRUE)")
    print("   • gene_dnds.parquet | bin_dnds.parquet")

if __name__ == "__main__":
    main()
