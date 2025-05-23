#!/usr/bin/env python3
# calculate_combined_stats_v8.py ------------------------------------------------
# This script processes mutation data to calculate various statistics and metrics
# for analyzing strain evolution over time. It handles:
# 1. Mutation frequency trajectories
# 2. Gene-level statistics
# 3. Functional enrichment (GO/KEGG terms)
# 4. dN/dS ratios
# 5. Sweeping variant detection
# ------------------------------------------------------------------------------

import argparse, re, sys, gzip, uuid, functools, operator as op
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds
from tqdm import tqdm

pd.options.mode.copy_on_write = True           # pandas ≥2.1

# ───────────────────────────── CLI ──────────────────────────────────────────
def get_args():
    """
    Parse command line arguments for the script.
    Required arguments:
    - mut-file: Input mutation data file
    - scaffold-file: File containing scaffold/bin mapping
    - meta-file: Metadata file with patient information
    Optional arguments:
    - out-dir: Output directory for results
    - chunk-size: Number of rows to process at once
    - low-mem: Use low memory mode
    - patient-chunks: Process data in patient-specific chunks
    - base-window: Baseline time window (e.g., "-2:0")
    - post-window: Post-treatment time window (e.g., "10:28")
    """
    p = argparse.ArgumentParser("iHMP mutation summariser (stream, sweeps, bins)")
    p.add_argument("--mut-file",      required=True)
    p.add_argument("--scaffold-file", required=True)
    p.add_argument("--meta-file",     required=True)
    p.add_argument("--out-dir",       default="summaries")
    p.add_argument("--chunk-size",    type=int, default=1_000_000)
    p.add_argument("--low-mem",       action="store_true")
    p.add_argument("--patient-chunks",action="store_true")
    p.add_argument("--base-window",   default="-2:0",
        help="baseline day range, e.g. -2:0 (inclusive)")
    p.add_argument("--post-window",   default="10:28",
        help="post-antibiotic day range, e.g. 10:28 (inclusive)")
    return p.parse_args()

def parse_win(win: str)->Tuple[int,int]:
    """
    Parse a window string (e.g., "-2:0") into start and end integers.
    Returns the minimum and maximum values to ensure correct ordering.
    """
    a,b = map(int, win.split(":"))
    return min(a,b), max(a,b)

# ────────────────────────── Helper Functions ─────────────────────────────────
def detect_week_cols(cols):
    """
    Identify numeric columns that represent weeks in the dataset.
    Returns a sorted list of week numbers.
    """
    return sorted(int(c) for c in cols if re.fullmatch(r"\d+", c))

def classify_dn_ds(mut_type):
    """
    Classify mutations as non-synonymous (dN) or synonymous (dS).
    Returns two boolean masks for dN and dS mutations.
    """
    m = mut_type.str.lower().fillna("")
    return m.isin(["missense","nonsense"]), m.eq("silent")


def split_dbxrefs(s):
    """
    Split database cross-references into a list of terms.
    Handles comma-separated values and strips whitespace.
    """
    return [t.strip().lower() for t in str(s).split(",") if t.strip()]

def extract_go(tok):
    """
    Extract GO (Gene Ontology) terms from a list of database references.
    """
    return [t for t in tok if t.startswith("go:")]

def extract_kegg(tok):
    """
    Extract KEGG terms from a list of database references.
    Includes both EC numbers and K numbers.
    """
    return [t for t in tok if t.startswith("ec:") or re.fullmatch(r"k\d{5}",t)]

# ───────────────────────── Bin Mapping ──────────────────────────────────────
def load_bin_map(scaffold_file, meta_file):
    """
    Create a mapping between patient-scaffold pairs and their corresponding bins.
    Combines information from scaffold file and metadata.
    Returns a dictionary with (patient_id, scaffold) tuples as keys and bin as value.
    """
    # Load and standardize metadata
    meta = pd.read_csv(meta_file, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(" ","_")
    sample_col = "external.id" if "external.id" in meta.columns else "external_id"
    meta = meta.rename(columns={sample_col:"sample"})[["sample","participant_id"]]\
             .rename(columns={"participant_id":"patient_id"}).dropna()

    scaf = pd.read_csv(scaffold_file, sep="\t",
                       usecols=["scaffold","Sample","bin"], low_memory=False)
    scaf.columns = scaf.columns.str.lower()
    mp = scaf.merge(meta, on="sample", how="left").dropna(subset=["patient_id"])
    return {(r.patient_id, r.scaffold): r.bin for r in mp.itertuples(index=False)}

# ───────────────────── Streaming Writer ────────────────────────────────────
def append_parquet(dir_path:Path, df):
    """
    Append a DataFrame to a Parquet dataset.
    Creates a new partition file for each chunk of data.
    """
    if df.empty: return
    dir_path.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   dir_path/f"part-{uuid.uuid4().hex}.parquet",
                   compression="snappy")

# ─────────── Incremental Aggregation ───────────────────────────────────────
def merge_stat(running, new, key_cols, agg_dict):
    """
    Merge new statistics with running totals.
    Performs incremental aggregation to handle large datasets efficiently.
    """
    if new.empty: return running
    combined = pd.concat([running,new]) if running is not None else new
    return combined.groupby(key_cols,as_index=False).agg(agg_dict)

# ─────────────────────── Chunk Processor ────────────────────────────────────
def process_chunk(df, week_cols, bin_map, base_rng, post_rng, control_ids:set):
    """
    Process a chunk of mutation data to calculate various statistics.
    
    Parameters:
    - df: Input DataFrame with mutation data
    - week_cols: List of week numbers
    - bin_map: Mapping of patient-scaffold pairs to bins
    - base_rng: Baseline time window
    - post_rng: Post-treatment time window
    - control_ids: Set of control patient IDs
    
    Returns:
    Tuple of DataFrames containing:
    1. Long-format frequency data
    2. SNP statistics
    3. Gene trajectory data
    4. Gene-level dN/dS ratios
    5. Bin-level dN/dS ratios
    6. GO term statistics
    7. KEGG term statistics
    8. Variant frequency ranges
    """
    # Convert numeric columns
    num_cols = ["p_value","slope","freq_range","min_freq","max_freq","r_sq",
                *map(str,week_cols)]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Add bin information
    if "bin" not in df.columns:
        df["bin"] = [bin_map.get((pid,chrom))
                     for pid,chrom in zip(df.patient_id,df.chromosome)]

    # Extract functional annotations
    if "dbxrefs" in df.columns and "go_terms" not in df.columns:
        toks = df.dbxrefs.apply(split_dbxrefs)
        df["go_terms"]   = toks.apply(extract_go)
        df["kegg_terms"] = toks.apply(extract_kegg)

    # Classify mutations
    mut_lower = df["mutation_type"].str.lower().fillna("")
    df["mut_class"] = np.select(
        [
            mut_lower.isin(["missense", "nonsense"]),
            mut_lower.eq("silent")
        ],
        ["coding", "silent"],
        default="intergenic"
    )

    # Filter out silent mutations for downstream analysis
    df_func = df[df.mut_class != "silent"].copy()
    df_ns  = df_func[df_func["mut_class"] == "coding"].copy()
    df_int = df_func[df_func["mut_class"] == "intergenic"].copy()

    # Convert to long format for frequency analysis
    id_cols = [c for c in df_func.columns if c not in map(str,week_cols)]
    long_df = (df_func.melt(id_vars=id_cols, value_vars=list(map(str,week_cols)),
                       var_name="week_num", value_name="freq").dropna(subset=["freq"]))
    long_df["week_num"] = long_df.week_num.astype(int)
    
    # Calculate per-variant frequency ranges
    var_range = (long_df.groupby(["patient_id", "chromosome", "position", "gene", "mutation_type"],
                             as_index=False)
                       .agg(min_f=("freq", "min"),
                            max_f=("freq", "max")))
 
    # Calculate SNP-level statistics
    snp_cols = [c for c in ["patient_id","chromosome","position",
                            "ref_base","new_base","slope","p_value",
                            "r_sq","min_freq","max_freq","freq_range"]
                if c in df_func.columns]
    snp_stats = df_func[snp_cols].drop_duplicates()

    def summarise_terms(df, group_cols):
        return (
            df.groupby(group_cols, as_index=False)
            .agg(n_snps=("p_value", "size"),
                min_p=("p_value", "min"),
                sum_p=("p_value", "sum"),
                sum_abs_slope=("slope", lambda x: np.nansum(np.abs(x))),
                sum_delta=("freq_range", "sum"))
        )

    # gene-level
    gene_ns  = summarise_terms(df_ns,  ["patient_id", "gene"]) if "gene" in df_ns else pd.DataFrame()
    gene_int = summarise_terms(df_int, ["patient_id", "gene"]) if "gene" in df_int else pd.DataFrame()

    # GO/KEGG-level
    go_ns = go_int = kegg_ns = kegg_int = pd.DataFrame()

    if "go_terms" in df_ns:
        go_ns = summarise_terms(df_ns.explode("go_terms").dropna(subset=["go_terms"]),
                                ["patient_id", "go_terms"])
    if "go_terms" in df_int:
        go_int = summarise_terms(df_int.explode("go_terms").dropna(subset=["go_terms"]),
                                ["patient_id", "go_terms"])
    if "kegg_terms" in df_ns:
        kegg_ns = summarise_terms(df_ns.explode("kegg_terms").dropna(subset=["kegg_terms"]),
                                ["patient_id", "kegg_terms"])
    if "kegg_terms" in df_int:
        kegg_int = summarise_terms(df_int.explode("kegg_terms").dropna(subset=["kegg_terms"]),
                                ["patient_id", "kegg_terms"])

    # Calculate dN/dS ratios
    gene_dnds = bin_dnds = pd.DataFrame()
    if "mutation_type" in df.columns:
        is_dn,is_ds = classify_dn_ds(df.mutation_type)
        # Gene-level dN/dS
        if "gene" in df.columns:
            gene_dnds = (df.assign(is_dn=is_dn,is_ds=is_ds)
                           .groupby(["patient_id","gene"],as_index=False)
                           .agg(dn=("is_dn","sum"), ds=("is_ds","sum")))
        # Bin-level dN/dS
        if "bin" in df.columns:
            bin_dnds = (df.assign(is_dn=is_dn,is_ds=is_ds)
                          .dropna(subset=["bin"])
                          .groupby(["patient_id","bin"],as_index=False)
                          .agg(dn=("is_dn","sum"), ds=("is_ds","sum")))

    return (long_df, snp_stats, gene_ns, gene_int, gene_dnds,
            bin_dnds, go_ns, go_int, kegg_ns, kegg_int, var_range)

# ───────────────────────────── Main ─────────────────────────────────────────
def main():
    """
    Main function that orchestrates the data processing pipeline:
    1. Parse command line arguments
    2. Load and prepare metadata
    3. Process data in chunks
    4. Calculate and aggregate statistics
    5. Write results to output files
    """
    # Parse arguments and setup
    a = get_args()
    out_dir = Path(a.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Parse time windows
    base_rng = parse_win(a.base_window)
    post_rng = parse_win(a.post_window)

    # Load metadata and identify control patients
    meta = pd.read_csv(a.meta_file, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(" ","_")
    control_ids = set(meta.loc[meta.antibiotics.isna(),"participant_id"])

    # Load bin mapping
    bin_map = load_bin_map(a.scaffold_file, a.meta_file)

    # Detect week columns in the data
    preview = pd.read_csv(a.mut_file, sep="\t", nrows=2000, low_memory=False)
    week_cols = detect_week_cols(preview.columns)
    if not week_cols:
        sys.exit("❌ No numeric week columns detected.")
    print(f"Detected {len(week_cols)} week columns; first 10 → {week_cols[:10]}")

    # Setup progress bar
    opener = gzip.open if a.mut_file.endswith(".gz") else open
    total_rows = sum(1 for _ in opener(a.mut_file,"rt")) - 1
    pbar = tqdm(total=total_rows, unit="rows")

    # Initialize running summary tables
    runs = {k: None for k in [
        "gene_ns", "gene_int", "gene_dnds", "bin_dnds",
        "go_ns", "go_int", "kegg_ns", "kegg_int",
        "sweep_core"
    ]}

    # Process data in chunks
    for chunk in pd.read_csv(a.mut_file, sep="\t", chunksize=a.chunk_size,
                             low_memory=not a.low_mem):
        # Standardize column names
        chunk.columns = (chunk.columns.str.strip()
                                      .str.lower()
                                      .str.replace(" ","_"))
        
        # Optionally split by patient
        sub_chunks = (chunk.groupby("participant_id",sort=False)
                            .apply(lambda g:g)) if a.patient_chunks else [chunk]
       
        for sub in chunk:
            if sub.empty: continue
            if "participant_id" in sub.columns:
                sub = sub.rename(columns={"participant_id":"patient_id"})
            if sub.patient_id.isna().all(): continue

            # Process chunk and get statistics
            dfs = process_chunk(sub, week_cols, bin_map,
                                base_rng, post_rng, control_ids)
            (long_df,snp_stats,gene_ns,gene_int,gene_dnds,
             bin_dnds,go_ns,go_int,kegg_ns,kegg_int,var_range) = dfs

            # Write large datasets to disk
            append_parquet(out_dir/"mutation_long", long_df)
            append_parquet(out_dir/"snp_stats",     snp_stats)

            # Update running statistics
            runs["sweep_core"] = merge_stat(
                    runs["sweep_core"], var_range,
                    ["patient_id", "chromosome", "position", "gene", "mutation_type"],
                    {"min_f": "min", "max_f": "max"})
            runs["gene_ns"] = merge_stat(runs["gene_ns"], gene_ns,
                ["patient_id","gene"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["gene_int"] = merge_stat(runs["gene_int"], gene_int,
                ["patient_id","gene"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["gene_dnds"]= merge_stat(runs["gene_dnds"], gene_dnds,
                ["patient_id","gene"],{"dn":"sum","ds":"sum"})
            runs["bin_dnds"] = merge_stat(runs["bin_dnds"], bin_dnds,
                ["patient_id","bin"],{"dn":"sum","ds":"sum"})
            runs["go_ns"] = merge_stat(runs["go_ns"], go_ns,
                ["patient_id","go_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["go_int"] = merge_stat(runs["go_int"], go_int,
                ["patient_id","go_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["kegg_ns"]=merge_stat(runs["kegg_ns"], kegg_ns,    
                ["patient_id","kegg_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["kegg_int"] = merge_stat(runs["kegg_int"], kegg_int,
                ["patient_id","kegg_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})

        pbar.update(len(chunk))
    pbar.close()

    # ── finalise universal sweep definition ──────────────────────────
    core = runs["sweep_core"]
    if core is not None:
        core["sweep"] = (core.min_f < 0.20) & (core.max_f > 0.80)
        pq.write_table(pa.Table.from_pandas(core, preserve_index=False),
                    out_dir / "sweep_stats.parquet",
                    compression="snappy")
        print(f"✔ sweep_stats.parquet ({core.shape[0]:,} variants)")


    # ── finalise other tables (weighted means, dN/dS) ───────────────────────
    def finalise_means(df,sum_cols):
        for col in sum_cols:
            df[col.replace("sum_","mean_")] = df[col]/df["n_snps"]
            df.drop(columns=[col],inplace=True)
        return df

    # after your existing "finalise_means()" loop ------------------------------
    def write_counts(df, key_cols, out_name):
        if df is None: return
        counts = df[key_cols + ["n_snps", "mean_abs_slope"]].copy()
        # Apply weighting: down-weight neutral mutations
        counts["weighted_n"] = counts["n_snps"] * (1 - counts["mean_abs_slope"].clip(upper=1)) # 1 - abs(slope)
        # Save to disk
        pq.write_table(pa.Table.from_pandas(counts[key_cols + ["weighted_n"]],
                                            preserve_index=False),
                    out_dir/f"{out_name}_counts.parquet",
                    compression="snappy")
        print(f"✔ {out_name}_counts.parquet ({counts.shape[0]:,} rows)")

    write_counts(runs["gene_traj"], ["patient_id","gene"],   "gene")
    write_counts(runs["kegg_stats"],["patient_id","kegg_terms"],"kegg")
    write_counts(runs["go_stats"],  ["patient_id","go_terms"], "go")

    for name,df in runs.items():
        if name=="sweep_core" or df is None: continue
        if name in ["gene_traj","go_stats","kegg_stats"]:
            df = finalise_means(df,["sum_abs_slope","sum_p","sum_delta"])
        if name.endswith("dnds"):
            df["dnds"] = df.apply(lambda r:r.dn/r.ds if r.ds else np.nan, axis=1)
        pq.write_table(pa.Table.from_pandas(df,preserve_index=False),
                       out_dir/f"{name}.parquet", compression="snappy")
        print(f"✔ {name}.parquet ({df.shape[0]:,} rows)")

    # ── burden per patient/week ─────────────────────────────────────────────
    long_ds = ds.dataset(out_dir/"mutation_long", format="parquet")
    long_df = long_ds.to_table(columns=["patient_id","week_num"]).to_pandas()
    burden  = (long_df.groupby(["patient_id","week_num"],as_index=False)
                        .size().rename(columns={"size":"n_snvs"}))
    pq.write_table(pa.Table.from_pandas(burden,preserve_index=False),
                   out_dir/"patient_week.parquet",compression="snappy")
    print(f"✔ patient_week.parquet ({burden.shape[0]:,} rows)")
    print("🏁 Finished — all summaries written to", out_dir)

if __name__=="__main__":
    main()
