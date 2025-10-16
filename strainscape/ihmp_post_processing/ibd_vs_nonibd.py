#!/usr/bin/env python3
"""
ibd_vs_nonibd.py

Streamlined IBD vs nonIBD analysis using Poisson GLM 

Model (per feature k):
  log E[Y_ik] = α + β * IBD_i + FE(MAG_i) + log(opportunity_ik)
  where Y_ik = # sweeping variants mapped to feature k in isolate i.

Key choices:
  • Opportunity (offset): ko_genes (default, auto-built) or non_sweep (fallback)
  • Fallback: Clustered Poisson GLM when FE saturates

Inputs
------
--parquet-dir        Path to parquet directory with SNV data
--meta-file          CSV/TSV with patient_id and diagnosis
--mag-genes-file     TSV with gene annotations
--filter-file        Optional CSV/TSV listing isolates to include
--opp-genes-file     Optional TSV with feature, n_genes

Main options
------------
--family             one of: gene (default), kegg, go, ec
--opp                opportunity type: ko_genes (default), non_sweep
--min-isolates       minimum # isolates per feature (default 10)
--sparse-cutoff      use binomial fallback if sweeps <= this (default 5)

Outputs
-------
out_dir/
  ibd_vs_nonibd_<family>.csv     # Feature-level results with RR, p, FDR, CIs
  summary.json                   # run stats

python strainscape/ihmp_post_processing/ibd_vs_nonibd.py \
    --parquet-dir /Users/reneeoles/Desktop/strainscape_output/iHMP/output/all_snvs \
    --meta-file /Users/reneeoles/Desktop/strainscape_output/iHMP/metadata/meta_map.csv \
    --mag-genes-file /Users/reneeoles/Desktop/strainscape_output/iHMP/input/all_genes.tsv \
    --filter-file /Users/reneeoles/Desktop/strainscape_output/iHMP/output/figures/sweeps/patient_bin.txt \
    --scaffold-file /Users/reneeoles/Desktop/strainscape_output/iHMP/input/combined_processed_scaffolds.txt \
    --out-dir /Users/reneeoles/Desktop/strainscape_output/iHMP/output/ibd_analysis_all \
    --family gene \
    --grouping dx3 \
    --min-isolates 10 \
    --background snv_blocks
"""

from __future__ import annotations
import argparse, json, logging, os, re, sys, time
from contextlib import contextmanager
from math import ceil
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm, t

try:
    import pyarrow.dataset as ds
    HAVE_ARROW = True
except Exception:
    HAVE_ARROW = False

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Poisson, Binomial
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.sandwich_covariance import cov_cluster
from statsmodels.tools.sm_exceptions import PerfectSeparationError

try:
    from joblib import Parallel, delayed
    HAVE_JOBLIB = True
except ImportError:
    HAVE_JOBLIB = False

# --- pandas dtype helpers (works with bool and nullable "boolean") ---
from pandas.api.types import is_bool_dtype

def _force_bool(s: pd.Series) -> pd.Series:
    """Return a plain bool dtype with NA→False, even if input is pandas 'boolean'."""
    # Convert pandas nullable boolean → fillna(False) → plain bool
    if str(s.dtype) == "boolean":
        return s.fillna(False).astype(bool)
    if is_bool_dtype(s):
        # already numpy bool_ dtype
        return s.astype(bool)
    # Fallback: anything truthy becomes True
    return s.astype(bool)


# --------------------------- logging ---------------------------------
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger("ibd_glm")

# --------------------------- timing helper ---------------------------
@contextmanager
def log_time(msg: str):
    t0 = time.perf_counter()
    log.info("%s …", msg)
    try:
        yield
    finally:
        log.info("%s ✓ (%.2fs)", msg, time.perf_counter() - t0)

# One delimiter regex for everything
DELIMS = re.compile(r"[,\s;|/]+")

# Compact, prefix-tolerant patterns
KO_RX = re.compile(r"^(?:kegg:)?(K\d{5})$", re.IGNORECASE)
GO_RX = re.compile(r"^(?:go:)?(\d{7})$", re.IGNORECASE)
EC_RX = re.compile(r"^(?:ec:)?(\d+(?:\.\d+){0,3})$", re.IGNORECASE)

GENE_BAD = {"hypothetical","hypothetical_protein","uncharacterized","putative","unknown","none","nan",""}

# Single-path constants (no configuration switches)
FREQ_CUT = 0.6      # frequency threshold for high-freq sweeps
MIN_EVENTS = int(os.environ.get("MIN_EVENTS_OVERRIDE", 10))  # total high-freq sweeps needed to fit
MIN_OPP_PER_GRP = 10  # per-group isolates with opportunity > 0
MIN_P = 1e-20       # minimum p-value floor for numerical stability


def _coerce_to_str(x) -> str:
    """Flatten np.ndarray/list to comma string; else str()."""
    if x is None or (isinstance(x, float) and np.isnan(x)): 
        return ""
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if isinstance(x, list):
        return ",".join(map(str, x))
    return str(x)

def two_sided_p_from_z(z: float, df: Optional[int] = None) -> float:
    """Convert z-score to two-sided p-value, optionally using t-distribution."""
    if df is None or df <= 0:
        # Normal distribution
        p = 2 * (1 - norm.cdf(abs(z)))
    else:
        # t-distribution with df degrees of freedom
        from scipy.stats import t
        p = 2 * (1 - t.cdf(abs(z), df))
    
    # Apply minimum p-value floor
    return max(p, MIN_P)

def has_minimum_effect_size(results: dict, min_effect_size: float = 0.5) -> bool:
    """Check if any group has minimum effect size (log2RR >= min_effect_size)."""
    for key, value in results.items():
        if key.startswith("RR_") and isinstance(value, (int, float)):
            log2rr = np.log2(value)
            if abs(log2rr) >= min_effect_size:
                return True
    return False

def _split_dbx(x):
    if pd.isna(x) or x is None: return []
    return [t for t in re.split(r"[,;|\s]+", str(x).strip()) if t and t.lower() not in {"nan","none"}]

def feats_from_dbxrefs(dbx, family):
    vals=[]
    for t in _split_dbx(dbx):
        u=t.upper()
        if family=="kegg":
            u = u.split(":",1)[1] if u.startswith("KEGG:") else u
            if KO_RX.match(u): vals.append(u)
        elif family=="go":
            u = u.split(":",1)[1] if u.startswith("GO:") else u
            if GO_RX.match(u): vals.append(u)
        elif family=="ec":
            u = u.split(":",1)[1] if u.startswith("EC:") else u
            if EC_RX.match(u): vals.append(u)
    return list(dict.fromkeys(vals))

def feats_from_gene(rep_gene):
    if pd.isna(rep_gene) or not str(rep_gene).strip(): return []
    toks = [t.strip() for t in re.split(r"[ ,;|]+", str(rep_gene)) if t.strip()]
    bad = {"hypothetical","hypothetical_protein","uncharacterized","putative","unknown","none","nan"}
    return [t for t in toks if t.lower() not in bad]

def parse_feature_list(x, family: str) -> List[str]:
    """Parse feature list using optimized functions."""
    if family == "gene":
        return feats_from_gene(x)
    else:
        return feats_from_dbxrefs(x, family)

def build_feature_map_from_mag(mag_genes: pd.DataFrame, family: str) -> pd.DataFrame:
    """Return columns: iso_id, locus_tag, feature (deduped)."""
    feats = (mag_genes["rep_gene"] if family=="gene"
             else mag_genes["rep_dbxrefs"]).apply(lambda x: parse_feature_list(x, family))
    fm = (mag_genes.assign(feature_list=feats)
            .explode("feature_list")
            .rename(columns={"feature_list":"feature"})
            .dropna(subset=["feature"])
            .loc[:, ["patient_id","bin","locus_tag","feature"]]
            .copy())
    fm["iso_id"] = fm["patient_id"].astype(str) + "|" + fm["bin"].astype(str)
    fm = fm.drop(columns=["patient_id","bin"]).drop_duplicates()
    return fm

def read_table_auto(p: Path) -> pd.DataFrame:
    """
    Truly minimal: let pandas infer delimiter.
    (Uses python engine; slightly slower but robust.)
    """
    return pd.read_csv(p, sep=None, engine="python")

def map_dx_to_group(meta: pd.DataFrame, grouping: str = "dx3") -> pd.DataFrame:
    """Normalize patient id + collapse diagnosis to nonIBD/UC/CD or nonIBD/IBD."""
    # Find/rename id column
    for c in ["patient_id","Participant ID","participant_id","Run","run","subject","sample"]:
            if c in meta.columns:
                meta = meta.rename(columns={c: "patient_id"})
                break
    if "patient_id" not in meta.columns or "diagnosis" not in meta.columns:
        raise ValueError("meta must have patient_id and diagnosis")
    
    dx = meta["diagnosis"].astype(str).str.strip()
    if grouping == "ibd2":
        mapping = {"CD":"IBD", "UC":"IBD", "nonIBD":"nonIBD"}
        cats = ["nonIBD","IBD"]
    else:
        mapping = {"CD":"CD", "UC":"UC", "nonIBD":"nonIBD"}
        cats = ["nonIBD","UC","CD"]

    out = (meta.assign(group=dx.map(mapping))
                 .dropna(subset=["group"])
                 .loc[:, ["patient_id","group"]]
                 .drop_duplicates("patient_id"))
    out["group"] = pd.Categorical(out["group"], categories=cats, ordered=False)
    return out

def build_opp_genes_from_mag(mag_genes: pd.DataFrame, family: str) -> pd.DataFrame:
    """
    Count features per MAG for opportunity (only if you keep this path).
    """
    feats = (mag_genes["rep_gene"] if family=="gene" else mag_genes["rep_dbxrefs"]).apply(
        lambda x: parse_feature_list(x, family)
    )
    fm = (mag_genes.assign(feature=feats)
                    .explode("feature")
                    .dropna(subset=["feature"]))
    opp = (fm.groupby(["feature"], observed=True)["locus_tag"]
             .nunique().rename("n_genes").reset_index())
    return opp

def read_snvs(parquet_dir: Path, cols_needed: List[str]) -> pd.DataFrame:
    """
    Fast parquet reader with PyArrow fallback to pandas.
    """
    with log_time("[LOAD] Reading SNVs"):
        if HAVE_ARROW:
            dataset = ds.dataset(str(parquet_dir))
            cols = [c for c in cols_needed if c != "*"]
            table = dataset.to_table(columns=cols)
            return table.to_pandas(types_mapper=pd.ArrowDtype)
    if parquet_dir.is_file():
        return pd.read_parquet(parquet_dir, columns=[c for c in cols_needed if c != "*"])
    parts = sorted(parquet_dir.glob("*.parquet"))
    if not parts:
        raise ValueError(f"No parquet files in {parquet_dir}")
    dfs = [pd.read_parquet(p, columns=[c for c in cols_needed if c != "*"]) for p in parts]
    return pd.concat(dfs, ignore_index=True)

def split_locus_tags(df: pd.DataFrame, col="locus_tag", weight_col="freq_range", sep="/") -> pd.DataFrame:
    """Optimized locus tag splitting from run_enrichment_analysis_clusters."""
    s = df[col].astype("string").str.split(sep)
    out = df.join(s.rename("_parts")).explode("_parts")
    out[col] = out["_parts"].astype("string").str.strip()
    out = out.drop(columns=["_parts"])
    out = out[out[col].notna() & (out[col]!="") & (~out[col].str.lower().isin({"nan","none"}))]
    if weight_col in out.columns:
        counts = s.str.len()
        out["_n"] = counts.loc[out.index].clip(lower=1)
        out[weight_col] = out[weight_col] / out["_n"]; out.drop(columns=["_n"], inplace=True)
    return out


def load_filter(filter_path: Optional[Path]) -> Optional[pd.DataFrame]:
    if not filter_path: return None
    flt = pd.read_csv(filter_path, sep=None, engine="python")
    req = ["patient_id","bin"]
    if not all(c in flt.columns for c in req):
        raise ValueError(f"Filter file needs columns: {req}")
    flt["patient_id"] = flt["patient_id"].astype(str).str.strip()
    flt["bin"] = flt["bin"].astype(str).str.strip()
    flt["iso_id"] = flt["patient_id"] + "|" + flt["bin"]
    return flt[["iso_id"]].drop_duplicates()

def apply_filter(df: pd.DataFrame, filter_file: Optional[Path]) -> pd.DataFrame:
    """Apply filter using optimized load_filter function."""
    if not filter_file:
        return df
    
    # Ensure iso_id exists on df
    if "iso_id" not in df.columns and {"patient_id","bin"}.issubset(df.columns):
        df = df.assign(iso_id=df["patient_id"].astype(str) + "|" + df["bin"].astype(str))
    
    filt = load_filter(filter_file)
    if filt is not None:
        keep = filt["iso_id"].unique()
        return df[df["iso_id"].isin(keep)]
    return df

def load_mag_genes(genes_file: Path) -> pd.DataFrame:
    """
    Load MAG genes from all_genes.tsv file
    
    Expected columns: patient_id, bin, scaffold, locus_tag, cluster_rep, rep_gene, rep_dbxrefs
    Returns: DataFrame with iso_id added
    """    
    # Load the genes data
    mag_genes = pd.read_csv(genes_file, sep='\t')
    
    # Create iso_id from patient_id and bin
    mag_genes["iso_id"] = mag_genes["patient_id"].astype(str) + "|" + mag_genes["bin"].astype(str)
    
    log.info("[LOAD] Unique iso_ids in MAG genes: %d", mag_genes["iso_id"].nunique())
    
    return mag_genes


def load_scaffold_coverage(scaffold_file: Path, meta_file: Path) -> pd.DataFrame:
    """
    Load scaffold coverage data from combined_processed_scaffolds.txt
    and map Run IDs to Participant IDs using metadata
    
    Expected columns: scaffold, length, coverage, breadth, Sample, bin, etc.
    Returns: iso_id, scaffold, breadth, coverage (with iso_id using Participant ID)
    """
    log.info("[LOAD] Loading scaffold coverage from: %s", scaffold_file)
    
    # Load the scaffold data
    scaffold_df = pd.read_csv(scaffold_file, sep='\t')
    
    # Load metadata to map Run -> Participant ID
    log.info("[LOAD] Loading metadata for Run -> Participant ID mapping from: %s", meta_file)
    meta = pd.read_csv(meta_file)
    run_to_participant = meta[["Run", "Participant ID"]].drop_duplicates()
    
    # Map Sample (Run) to Participant ID
    scaffold_df = scaffold_df.merge(run_to_participant, left_on="Sample", right_on="Run", how="left")
    
    # Create iso_id from Participant ID and bin (to match MAG genes format)
    scaffold_df["iso_id"] = scaffold_df["Participant ID"].astype(str) + "|" + scaffold_df["bin"].astype(str)
    
    # Drop rows where Participant ID mapping failed
    scaffold_df = scaffold_df.dropna(subset=["Participant ID"])
    
    # Select and rename relevant columns
    scaffold_cov = scaffold_df[["iso_id", "scaffold", "breadth", "coverage"]].copy()
    
    log.info("[LOAD] Loaded scaffold coverage: %d scaffolds, %d unique isolates", 
             len(scaffold_cov), scaffold_cov["iso_id"].nunique())
    log.info("[LOAD] Coverage range: %.2f - %.2f, Breadth range: %.3f - %.3f", 
             scaffold_cov["coverage"].min(), scaffold_cov["coverage"].max(),
             scaffold_cov["breadth"].min(), scaffold_cov["breadth"].max())
    
    return scaffold_cov

def extract_gene_lengths_from_snvs(snvs: pd.DataFrame,
                                   valid_locus_tags: Optional[pd.Series] = None) -> pd.DataFrame:
    # subset to needed cols
    cols = ["locus_tag", "start", "stop"]
    sn = snvs.loc[:, [c for c in cols if c in snvs.columns]].dropna(subset=["locus_tag", "start", "stop"]).copy()
    if valid_locus_tags is not None:
        sn = sn[sn["locus_tag"].isin(valid_locus_tags)]

    s_start = sn["start"].astype(str)
    s_stop  = sn["stop"].astype(str)

    both_slash = s_start.str.contains("/", regex=False, na=False) & s_stop.str.contains("/", regex=False, na=False)
    neither    = ~s_start.str.contains("/", regex=False, na=False) & ~s_stop.str.contains("/", regex=False, na=False)
    mixed      = ~(neither | both_slash)

    L = pd.Series(np.nan, index=sn.index, dtype="float64")

    # simple: stop - start + 1
    if neither.any():
        n_start = pd.to_numeric(s_start[neither], errors="coerce")
        n_stop  = pd.to_numeric(s_stop[neither],  errors="coerce")
        simple = n_stop - n_start + 1
        L.loc[neither] = simple.to_numpy()

    # intergenic: start1 - stop0 + 1
    if both_slash.any():
        s_parts = s_start[both_slash].str.split("/", n=1, expand=True)
        e_parts = s_stop[both_slash].str.split("/", n=1, expand=True)
        s1 = pd.to_numeric(s_parts[1], errors="coerce")
        e0 = pd.to_numeric(e_parts[0], errors="coerce")
        inter = s1 - e0 + 1
        L.loc[both_slash] = inter.to_numpy()

    # mixed fallback: use first tokens
    if mixed.any():
        ms = s_start[mixed].str.split("/", n=1, expand=True)[0]
        me = s_stop[mixed].str.split("/", n=1, expand=True)[0]
        m_start = pd.to_numeric(ms, errors="coerce")
        m_stop  = pd.to_numeric(me, errors="coerce")
        mlen = m_stop - m_start + 1
        L.loc[mixed] = mlen.to_numpy()

    # clean invalid values
    L = L.where(np.isfinite(L) & (L > 0))

    out = (sn.assign(gene_len=L)
             .dropna(subset=["gene_len"])
             .groupby("locus_tag", observed=True)["gene_len"].median()
             .reset_index())
    return out

def make_callable_genes(
    mag_genes: pd.DataFrame,
    scaffold_cov: pd.DataFrame,
    snvs: pd.DataFrame,
    min_breadth: float = 0.80,
    min_depth: float = 10.0,
    default_gene_len: int = 1000,
    use_snv_lens: bool = True
) -> pd.DataFrame:
    """
    Build per-iso×gene callability from scaffold-level coverage/breadth.

    mag_genes: columns must include patient_id, bin, scaffold, locus_tag
    scaffold_cov: columns iso_id, scaffold, breadth (0-1), coverage (depth)
    snvs: SNV data to extract gene lengths from
    default_gene_len: fallback length for genes not found in SNV data

    Returns: iso_id, locus_tag, is_callable, callable_bases, gene_len
    """
    mg = mag_genes.copy()
    mg["iso_id"] = mg["patient_id"].astype(str) + "|" + mg["bin"].astype(str)

    # Extract gene lengths from SNV data (vectorized, fast) or use default
    valid_tags = mag_genes["locus_tag"].dropna().unique()
    if use_snv_lens:
        log.info("[CALLABLE] Extracting gene lengths from SNV data...")
        gene_lengths = extract_gene_lengths_from_snvs(snvs, valid_locus_tags=valid_tags)
        log.info("[CALLABLE] Extracted lengths for %d genes from SNV data", len(gene_lengths))
    else:
        log.info("[CALLABLE] Using default gene length (%d bp) for all genes", default_gene_len)
        gene_lengths = pd.DataFrame({"locus_tag": valid_tags, "gene_len": default_gene_len})
    
    # Merge gene lengths with mag_genes
    mg = mg.merge(gene_lengths, on="locus_tag", how="left")
    
    # Fill missing lengths with default
    mg["gene_len"] = mg["gene_len"].fillna(default_gene_len)
    
    log.info("[CALLABLE] Gene length statistics: mean=%.1f, median=%.1f, range=%.0f-%.0f", 
             mg["gene_len"].mean(), mg["gene_len"].median(), 
             mg["gene_len"].min(), mg["gene_len"].max())

    cov = scaffold_cov[["iso_id","scaffold","breadth","coverage"]].copy()
    cc = mg.merge(cov, on=["iso_id","scaffold"], how="left")

    cc["is_callable"] = (cc["breadth"] >= float(min_breadth)) & (cc["coverage"] >= float(min_depth))
    cc["is_callable"] = _force_bool(cc["is_callable"])
    cc["callable_bases"] = np.where(
        cc["is_callable"],
        cc["gene_len"].astype(float) * cc["breadth"].astype(float),
        0.0
    )
    
    log.info("[CALLABLE] Created callable genes: %d gene-isolate pairs, %d callable (%.1f%%)", 
             len(cc), cc["is_callable"].sum(),
             100 * cc["is_callable"].sum() / len(cc) if len(cc) > 0 else 0)
    
    return cc[["iso_id","locus_tag","is_callable","callable_bases","gene_len"]]

def _detect_pos_col(df: pd.DataFrame) -> str | None:
    for c in ("pos","position","site","start","coord"):
        if c in df.columns:
            return c
    return None

# Position parsing regex for extracting first numeric token
POS_NUM_RX = re.compile(r"([0-9]+(?:\.[0-9]+)?)")

def _parse_pos_series(s: pd.Series) -> pd.Series:
    """Parse position series, taking the first numeric token (works for '150723/152009' etc.)."""
    return pd.to_numeric(
        s.astype(str).str.extract(POS_NUM_RX, expand=False),
        errors="coerce"
    )

def _parse_position(pos_str):
    """Parse position string, handling intergenic regions like '150723/152009'."""
    try:
        if pd.isna(pos_str):
            return np.nan
        
        pos_str = str(pos_str)
        
        # Simple position (e.g., "150723")
        if "/" not in pos_str:
            return int(pos_str)
        
        # Intergenic position (e.g., "150723/152009")
        # Use the midpoint for clustering purposes
        parts = pos_str.split("/")
        if len(parts) >= 2:
            start = int(parts[0])
            end = int(parts[1])
            return int((start + end) / 2)
        else:
            return int(parts[0])
            
    except (ValueError, IndexError):
        return np.nan


def assign_block_ids(sn: pd.DataFrame, gcols: List[str], pos_col: Optional[str],
                     max_bp_gap: int, max_freq_diff: float) -> pd.Series:
    """
    Vectorized single-linkage cluster IDs per (iso_id,locus_tag) sorted by position.
    A new block starts when:
      - group changes, OR
      - pos gap > max_bp_gap, OR
      - |Δfreq| > max_freq_diff
    Rows without pos get block 0.
    """
    if pos_col is None or pos_col not in sn.columns:
        return pd.Series(0, index=sn.index, dtype="Int64")

    d = sn.dropna(subset=[pos_col]).copy()
    if d.empty:
        return pd.Series(0, index=sn.index, dtype="Int64")

    # Sort and compute diffs
    d.sort_values(gcols + [pos_col], inplace=True)
    grp_change = ~d[gcols].eq(d[gcols].shift()).all(axis=1).fillna(True)
    pos  = d[pos_col].to_numpy()
    freq = d["freq_range"].to_numpy()
    pos_diff  = np.r_[np.inf, np.diff(pos)]
    freq_diff = np.r_[np.inf, np.abs(np.diff(freq))]

    new_block = grp_change | (pos_diff > max_bp_gap) | (freq_diff > max_freq_diff)
    block_ids = np.cumsum(new_block) - 1  # 0-based per entire sorted d

    # Rebase block ids to be contiguous within each (iso,locus_tag)
    # by subtracting the starting id for each group.
    grp_keys = d[gcols].astype(str).agg("|".join, axis=1)
    grp_first = d.groupby(grp_keys)[d.columns[0]].apply(lambda x: x.index[0])
    start_id = pd.Series(block_ids, index=d.index).groupby(grp_keys).transform("min")
    block_ids_local = block_ids - start_id.to_numpy()

    out = pd.Series(0, index=sn.index, dtype="Int64")
    out.loc[d.index] = block_ids_local
    return out


# Start algorithms 
def build_agg_snv_blocks(
    snvs: pd.DataFrame,
    meta: pd.DataFrame,
    family: str,
    mag_genes: pd.DataFrame,
    callable_genes: pd.DataFrame | None = None,
    # block / LD parameters
    freq_cut: float = 0.6,
    block_size_bp: int = 200,
    max_bp_gap: int | None = None,
    max_freq_diff: float = 0.15,
) -> pd.DataFrame:
    """
    Collapse SNVs to LD-like blocks per iso×gene, then aggregate to iso×feature.
    Denominator n = callable blocks (preferred: ceil(callable_len / block_size_bp);
    fallback: 1 per callable gene; if callable_genes is None, assume all genes callable).
    """
    sn = snvs.copy()
    sn["iso_id"] = sn["patient_id"].astype(str) + "|" + sn["bin"].astype(str)
    sn["is_high_site"] = (sn["freq_range"] >= float(freq_cut))

    pos_col = _detect_pos_col(sn)
    if pos_col is not None:
        sn["pos_bp"] = _parse_pos_series(sn[pos_col])
        if "stop" in sn.columns:
            # fill gaps from stop when start couldn't be parsed
            sn["pos_bp"] = sn["pos_bp"].fillna(_parse_pos_series(sn["stop"]))
        pos_col = "pos_bp"          # use the clean numeric column
        log.info("[SNV_BLOCKS] Parsed positions from %s, %d valid positions", pos_col, sn["pos_bp"].notna().sum())
    else:
        pos_col = None              # no coords found → one block per gene
        log.info("[SNV_BLOCKS] No position column found, will create one block per gene")
    
    max_bp_gap = int(max_bp_gap if max_bp_gap is not None else block_size_bp)

    # --- Feature map (gene -> feature), same parsing you already use ---
    fm = build_feature_map_from_mag(mag_genes, family)

    # --- Optional: callable genes with optional callable_len ---
    if callable_genes is None:
        callable_genes = (
            mag_genes[["patient_id","bin","locus_tag"]]
            .assign(iso_id=lambda x: x["patient_id"].astype(str)+"|"+x["bin"].astype(str),
                    is_callable=True,
                    callable_bases=np.nan)  # placeholder
            [["iso_id","locus_tag","is_callable","callable_bases"]]
        )

    # ensure columns exist
    if "is_callable" not in callable_genes.columns:
        callable_genes["is_callable"] = True
    if "callable_bases" not in callable_genes.columns:
        callable_genes["callable_bases"] = np.nan

    # --- Block assignment within iso×gene ---
    gcols = ["iso_id","locus_tag"]
    with log_time("[SNV_BLOCKS] Block assignment"):
        sn["_block_id"] = assign_block_ids(sn, gcols, pos_col, max_bp_gap, max_freq_diff)
    log.info("[SNV_BLOCKS] Unique blocks: %d (%.1f%% rows with NA pos -> block 0)",
             sn["_block_id"].nunique(dropna=True),
             100.0 * (sn[pos_col].isna().mean() if pos_col else 0.0))

    # --- Summarize blocks per iso×gene ---
    # Any site in block high → block is high
    block_hi = (sn.groupby(gcols + ["_block_id"], observed=True)["is_high_site"]
                    .any()
                  .rename("block_is_high")
                    .reset_index())

    blocks_per_gene = (block_hi.groupby(gcols, observed=True)
                              .agg(n_blocks_obs=("block_is_high","size"),
                                   y_blocks_hi=("block_is_high","sum"))
             .reset_index())

    # merge callable; estimate callable blocks
    log.info("[SNV_BLOCKS] Merging callable genes info with %d block summaries", len(blocks_per_gene))
    b = (blocks_per_gene
         .merge(callable_genes[["iso_id","locus_tag","is_callable","callable_bases"]],
                on=gcols, how="left")
         .rename(columns={"callable_bases":"callable_len"}))
    
    log.info("[SNV_BLOCKS] After merge: %d iso×gene combinations, %d callable (%.1f%%)", 
             len(b), b["is_callable"].sum(), 100 * b["is_callable"].mean())

    # n_callable_blocks:
    #   preferred: ceil(callable_len / block_size_bp) if callable & len known
    #   fallback: 1 if callable else 0
    # Vectorized calculation of n_callable_blocks
    cl = b["callable_len"].astype(float)
    is_call = b["is_callable"].fillna(False).astype(bool)
    n_blocks = np.where(is_call & cl.gt(0),
                        np.ceil(cl / float(block_size_bp)).clip(lower=1),
                        np.where(is_call, 1, 0))
    b["n_callable_blocks"] = n_blocks.astype(int)
    # cap y by n just in case
    b["y_blocks_hi"] = np.minimum(b["y_blocks_hi"].astype(int), b["n_callable_blocks"])
    
    log.info("[SNV_BLOCKS] Callable blocks calculated: mean=%.1f, median=%.1f, range=%d-%d", 
             b["n_callable_blocks"].mean(), b["n_callable_blocks"].median(),
             b["n_callable_blocks"].min(), b["n_callable_blocks"].max())

    # --- Map to features (feature-specific opportunities) ---
    log.info("[SNV_BLOCKS] Mapping blocks to features...")
    bf = (b.merge(fm, on=gcols, how="inner")
            .loc[:, ["iso_id","locus_tag","feature","y_blocks_hi","n_callable_blocks"]])
    
    log.info("[SNV_BLOCKS] After feature mapping: %d iso×gene×feature combinations", len(bf))

    agg = (bf.groupby(["iso_id","feature"], observed=True)
            .agg(y=("y_blocks_hi","sum"),
                 n=("n_callable_blocks","sum"))
            .reset_index())
    
    log.info("[SNV_BLOCKS] Final aggregation: %d iso×feature combinations", len(agg))

    # --- Exposure E: total high blocks across all features per isolate ---
    iso_E = (b.groupby("iso_id", observed=True)["y_blocks_hi"]
               .sum().rename("E").astype(int).reset_index())
    agg = agg.merge(iso_E, on="iso_id", how="left")
    
    # Create E_excl to avoid over-adjustment (exclude current feature's contribution)
    agg["E_excl"] = (agg["E"] - agg["y"]).clip(lower=0)
    log.info("[SNV_BLOCKS] Exposure statistics: E mean=%.1f, E_excl mean=%.1f", 
             agg["E"].mean(), agg["E_excl"].mean())

    # attach group
    pid_to_group = (meta.assign(patient_id=meta["patient_id"].astype(str))
                         .set_index("patient_id")["group"])
    agg["group"] = agg["iso_id"].astype(str).str.split("|", n=1).str[0].map(pid_to_group)

    # validity
    agg = agg[(agg["n"] > 0)]
    bad = agg[agg["y"] > agg["n"]]
    if not bad.empty:
        # very conservative fix (should be rare): clip
        agg.loc[agg["y"] > agg["n"], "y"] = agg.loc[agg["y"] > agg["n"], "n"]

    agg["y"] = agg["y"].astype(int)
    agg["n"] = agg["n"].astype(int)
    agg["E"] = agg["E"].fillna(0).astype(int)
    
    # Create E_excl to avoid over-adjustment (exclude current feature's contribution)
    agg["E_excl"] = (agg["E"] - agg["y"]).clip(lower=0)
    log.info("[GENE] Exposure statistics: E mean=%.1f, E_excl mean=%.1f", 
             agg["E"].mean(), agg["E_excl"].mean())

    # Convert to categoricals for memory efficiency
    for c in ("feature","iso_id"):
        if c in agg.columns:
            agg[c] = agg[c].astype("category")
    if "group" in agg.columns and not pd.api.types.is_categorical_dtype(agg["group"]):
        agg["group"] = pd.Categorical(agg["group"], categories=["nonIBD","UC","CD"], ordered=False)

    return agg[["iso_id","group","feature","y","n","E","E_excl"]]


def build_agg_gene(
    snvs: pd.DataFrame,
    meta: pd.DataFrame,
    family: str,
    mag_genes: pd.DataFrame,
    freq_cut: float = 0.6,
    callable_genes: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Output: iso_id, group, feature, y, n, E
    y = # swept callable genes in feature k
    n = # callable genes annotated to feature k
    E = total # swept callable genes across all features in the isolate
    """
    sn = snvs.copy()
    sn["iso_id"] = sn["patient_id"].astype(str) + "|" + sn["bin"].astype(str)
    sn["is_high_site"] = (sn["freq_range"] >= float(freq_cut))

    # Any swept gene (within callable set later)
    swept_gene = (sn.groupby(["iso_id","locus_tag"], observed=True)["is_high_site"]
                    .any().rename("is_swept_gene").reset_index())

    # Parse features from MAG annotations
    fm = build_feature_map_from_mag(mag_genes, family)
    # Add scaffold column if present in mag_genes
    if "scaffold" in mag_genes.columns:
        scaffold_map = mag_genes[["iso_id", "locus_tag", "scaffold"]].drop_duplicates()
        fm = fm.merge(scaffold_map, on=["iso_id", "locus_tag"], how="left")

    # ---- NEW: restrict to callable genes ----
    if callable_genes is not None:
        req = {"iso_id","locus_tag","is_callable"}
        if not req.issubset(callable_genes.columns):
            raise ValueError(f"callable_genes missing required columns {req}")

        cg = callable_genes[["iso_id","locus_tag","is_callable"]].copy()
        # normalize to plain bool (handles pandas 'boolean' as well)
        cg["is_callable"] = _force_bool(cg["is_callable"])

        fm = (fm.merge(cg, on=["iso_id","locus_tag"], how="left")
                .query("is_callable == True")
                .drop(columns=["is_callable"]))

        swept_gene = (swept_gene.merge(cg, on=["iso_id","locus_tag"], how="left")
                                .query("is_callable == True")
                                .drop(columns=["is_callable"]))
    else:
        log.warning("Proceeding without callable-gene filter — denominator inflation likely.")

    # Trials n: callable genes in feature
    n_tbl = (fm.groupby(["iso_id","feature"], observed=True)["locus_tag"]
               .nunique().rename("n").reset_index().astype({"n": int}))

    # Successes y: swept callable genes in feature
    merged = swept_gene.merge(
        fm[["iso_id","locus_tag","feature"]],
        on=["iso_id","locus_tag"], how="inner"
    )

    # Force boolean dtype; anything truthy/1 stays True, NaNs->False
    merged["is_swept_gene"] = merged["is_swept_gene"].astype(bool)

    y_tbl = (merged.loc[merged["is_swept_gene"]]
                   .groupby(["iso_id","feature"], observed=True)["locus_tag"]
             .nunique()
             .rename("y")
             .reset_index()
             .astype({"y": int}))

    # Exposure E: total swept callable genes per isolate
    iso_E = (swept_gene["is_swept_gene"]
             .groupby(swept_gene["iso_id"], observed=True)
             .sum().rename("E").astype(int).reset_index())

    agg = (n_tbl.merge(y_tbl, on=["iso_id","feature"], how="left")
                 .merge(iso_E, on="iso_id", how="left"))
    agg["y"] = agg["y"].fillna(0).astype(int)
    agg["E"] = agg["E"].fillna(0).astype(int)
    
    # Create E_excl to avoid over-adjustment (exclude current feature's contribution)
    agg["E_excl"] = (agg["E"] - agg["y"]).clip(lower=0)

    # Attach group
    pid_to_group = (meta.assign(patient_id=meta["patient_id"].astype(str))
                         .set_index("patient_id")["group"])
    agg["group"] = agg["iso_id"].str.split("|", n=1).str[0].map(pid_to_group)

    # Validity
    agg = agg[(agg["n"] > 0) & (agg["y"] >= 0) & (agg["y"] <= agg["n"])].copy()
    
    # Convert to categoricals for memory efficiency
    for c in ("feature","iso_id"):
        if c in agg.columns:
            agg[c] = agg[c].astype("category")
    if "group" in agg.columns and not pd.api.types.is_categorical_dtype(agg["group"]):
        agg["group"] = pd.Categorical(agg["group"], categories=["nonIBD","UC","CD"], ordered=False)
    
    return agg[["iso_id","group","feature","y","n","E","E_excl"]]


def _winsorize(x, p=0.995):
    if len(x) == 0:
        return x
    hi = np.quantile(x, p)
    return np.clip(x, 0, hi)

def fit_one_feature(dfk: pd.DataFrame, alpha: float = 1e-3, use_E_excl: bool = True, no_E: bool = False) -> Optional[dict]:
    if dfk.empty:
        return None
    d = dfk.copy()

    # Valid rows
    d = d[(d["n"] > 0) & (d["y"] >= 0) & (d["y"] <= d["n"])].copy()
    if d["y"].sum() < MIN_EVENTS:
        return None

    # Categorical group with nonIBD as reference
    d["group"] = pd.Categorical(d["group"].astype(str), categories=["nonIBD","UC","CD"], ordered=False)

    # Sanity guards before modeling
    d = d[d["group"].notna() & d["n"].gt(0)]
    if d["group"].nunique() < 2:
        return None

    # --- Jeffreys smoothing: mitigates separation at row level ---
    d["y_s"] = d["y"].astype(float) + 0.5
    d["n_s"] = d["n"].astype(float) + 1.0
    # freq_weights = n_s gives each row the right binomial weight
    fw = d["n_s"].astype(float)
    if not np.all(np.isfinite(fw)):
        fw = None
    else:
        # Cap weights to reduce leverage from very large n_s
        fw = np.minimum(fw, np.quantile(fw, 0.99))

    # --- E handling: winsorize + center log1p(E) ---
    if no_E:
        formula = "I(y_s / n_s) ~ C(group, Treatment(reference='nonIBD'))"
        E_col = "none"
    else:
        E_col = "E_excl" if (use_E_excl and "E_excl" in d.columns) else "E"
        E = d[E_col].astype(float).clip(lower=0)
        Ew = _winsorize(E.values, p=0.995)
        logE_c = np.log1p(Ew)
        logE_c = logE_c - np.nanmean(logE_c)
        d["logE_c"] = logE_c
        formula = "I(y_s / n_s) ~ C(group, Treatment(reference='nonIBD')) + logE_c"

    # Binomial GLM on smoothed proportion with weights = n_s
    mod = smf.glm(
        formula=formula,
        data=d,
        family=sm.families.Binomial(),
        freq_weights=fw
    )

    # Fit with quasi-binomial for overdispersion; if unstable, escalate ridge (L2) until sane
    try:
        res = mod.fit(method="newton", maxiter=100, disp=0, scale="X2")  # quasi-binomial
        params, bse = res.params, res.bse
        if np.any(~np.isfinite(bse)) or np.any(np.abs(params) > 15):
            raise RuntimeError("unstable MLE")
    except Exception:
        lam = alpha
        for _ in range(6):
            res_pen = mod.fit_regularized(alpha=lam, L1_wt=0.0, refit=False)
            params = res_pen.params
            if np.all(np.isfinite(params)) and np.all(np.abs(params) <= 12):
                # Try to refit from penalized start to recover SEs
                try:
                    res = mod.fit(start_params=params, method="newton", maxiter=100, disp=0)
                    params, bse = res.params, res.bse
                except Exception:
                    res = res_pen
                    bse = np.full_like(params, np.nan, dtype=float)
                break
            lam *= 5.0
        else:
            # still unstable → give up on this feature
            return None

    # Cluster-robust SEs by patient
    groups = d["iso_id"].astype(str).str.split("|", n=1).str[0]
    try:
        rob = res.get_robustcov_results(cov_type="cluster", groups=groups)
        params, bse = rob.params, rob.bse
        phi = float(rob.scale)  # quasi-binomial dispersion parameter
    except Exception:
        params, bse = res.params, getattr(res, "bse", res.bse)
        phi = float(getattr(res, "scale", 1.0))  # fallback to model scale

    out = {
        "effect": "Binomial_Jeffreys_logE_ridge",
        "phi": phi,  # quasi-binomial dispersion parameter
        "E_used": E_col,  # which exposure column was used
        "n_isolates_used": int(d["iso_id"].nunique()),
        "n_patients_used": int(groups.nunique()),
        "rows": int(len(d)),
        "mean_n": float(d["n"].mean()),
        "mean_E": float(d["E"].mean()),
        "sum_y_nonIBD": float(d.loc[d["group"]=="nonIBD","y"].sum()),
        "sum_n_nonIBD": float(d.loc[d["group"]=="nonIBD","n"].sum()),
        "sum_E_nonIBD": float(d.loc[d["group"]=="nonIBD","E"].sum()),
        "sum_y_UC": float(d.loc[d["group"]=="UC","y"].sum()),
        "sum_n_UC": float(d.loc[d["group"]=="UC","n"].sum()),
        "sum_E_UC": float(d.loc[d["group"]=="UC","E"].sum()),
        "sum_y_CD": float(d.loc[d["group"]=="CD","y"].sum()),
        "sum_n_CD": float(d.loc[d["group"]=="CD","n"].sum()),
        "sum_E_CD": float(d.loc[d["group"]=="CD","E"].sum()),
    }

    # Extract contrasts vs nonIBD
    for lev in ["UC","CD"]:
        term = f"C(group, Treatment(reference='nonIBD'))[T.{lev}]"
        if term in params.index:
            beta = float(params[term])
            se   = float(bse[params.index.get_loc(term)]) if hasattr(bse, "__len__") else np.nan
            z    = beta / se if (np.isfinite(se) and se > 0) else np.nan
            p    = two_sided_p_from_z(z)
            out[f"OR_{lev}"] = float(np.exp(beta)) if np.isfinite(beta) else np.nan
            out[f"p_{lev}"]  = p

    return out



# --------------------------- main (refactored) ---------------------------------
def _fit_one_named(feat, dko, use_E_excl, no_E):
    """Helper for parallel fitting."""
    info = fit_one_feature(dko, use_E_excl=use_E_excl, no_E=no_E)
    if info: 
        info["feature"] = feat
    return info

def _fit_all_features(agg: pd.DataFrame, min_isolates: int, use_E_excl: bool = True, no_E: bool = False) -> pd.DataFrame:
    """Filter features to those with enough isolates and fit per-feature models."""
    log.info("[FILTER] Starting feature filtering with min_isolates=%d", min_isolates)
    
    # Diagnostics
    diag = agg.groupby("feature")["iso_id"].nunique().rename("n_obs").reset_index()
    gcover = (agg.groupby(["feature","group"])["iso_id"].nunique()
                .unstack(fill_value=0).reset_index())
    
    log.info("[FILTER] Initial feature diagnostics: %d features total", len(diag))

    # Base thresholds
    diag2 = (diag.loc[(diag["n_obs"] >= min_isolates)]
                  .merge(gcover, on="feature", how="left"))
    
    log.info("[FILTER] After min_isolates filter (≥%d): %d features remain", min_isolates, len(diag2))

    # Each group must occur in ≥ 1 MAG (relaxed for small test datasets)
    static_cols = {"feature","n_obs"}
    group_cols  = [c for c in diag2.columns if c not in static_cols]
    diag2["ok_groups"] = (diag2[group_cols] >= 1).all(axis=1) if group_cols else False
    keep_feats = set(diag2.loc[diag2["ok_groups"], "feature"])

    # Filter aggregation
    agg_before = len(agg)
    agg = agg[agg["feature"].isin(keep_feats)].copy()
    log.info("[FILTER] Aggregation after feature filtering: %d -> %d rows (%.1f%% retained)", 
             agg_before, len(agg), 100 * len(agg) / agg_before if agg_before > 0 else 0)
    
    if not len(keep_feats):
        raise SystemExit("No features left after filtering; relax thresholds.")
    agg = agg[agg["feature"].isin(keep_feats)].copy()

    events = agg.groupby("feature", observed=True)["y"].sum().rename("total_y").reset_index()
    agg = agg.merge(events, on="feature", how="left")
    agg = agg[agg["total_y"] >= MIN_EVENTS].drop(columns=["total_y"])
    if agg.empty:
        raise SystemExit("No features have ≥ MIN_EVENTS sweeping genes after filters.")


    # Fit
    features_to_fit = agg["feature"].nunique()
    pairs = list(agg.groupby("feature", sort=False))
    
    if getattr(argparse, "_parsed_args", None) and getattr(argparse._parsed_args, "n_jobs", 1) > 1 and HAVE_JOBLIB:
        n_jobs = argparse._parsed_args.n_jobs
        with log_time(f"[FIT] Parallel fitting ({n_jobs} jobs)"):
            results = Parallel(n_jobs=n_jobs, prefer="processes")(
                delayed(_fit_one_named)(feat, dko, use_E_excl, no_E) for feat, dko in pairs
            )
        results = [r for r in results if r is not None]
    else:
        # Serial fallback
        results = []
        with log_time(f"[FIT] Fitting {features_to_fit} features"):
            for i, (feat, dko) in enumerate(pairs, 1):
                if i % 100 == 0 or i == features_to_fit:
                    log.info("[FIT] %d/%d: %s", i, features_to_fit, feat)
                info = fit_one_feature(dko, use_E_excl=use_E_excl, no_E=no_E)
                if info is None:
                    continue
                info["feature"] = feat
                results.append(info)
    
    if not results:
        raise SystemExit("No features produced a valid fit. Try lowering filters.")
    res = pd.DataFrame(results)
    return res, agg

def _add_group_cc_pseudorows(d: pd.DataFrame, wt=0.5):
    # d has columns: y, n, E, group, iso_id
    sums = d.groupby("group", observed=True).agg(y=("y","sum"), n=("n","sum"))
    need_cc = (sums["y"]==0) | (sums["y"]==sums["n"])
    if not need_cc.any():
        return d

    medE = float(np.median(d["E"].clip(lower=1)))
    rows = []
    for g, flag in need_cc.items():
        if not flag: 
            continue
        rows.append({
            "iso_id": f"pseudo|{g}",
            "group": g,
            "y": 0.5,          # Haldane–Anscombe
            "n": 1.0,
            "E": medE,
            "offset_logE": np.log(medE),
            "wt": wt
        })
    if not rows:
        return d
    d2 = pd.concat([d, pd.DataFrame(rows)], ignore_index=True)
    # ensure categories unchanged
    if pd.api.types.is_categorical_dtype(d["group"]):
        d2["group"] = d2["group"].astype(d["group"].dtype)
    return d2


def _compute_fdrs_and_effects(res: pd.DataFrame) -> pd.DataFrame:
    res = res.copy()
    p_cols = [c for c in res.columns if c.startswith("p_")]
    groups_seen = set()

    if p_cols:
        for pc in p_cols:
            grp = pc.replace("p_", "")
            p = pd.to_numeric(res[pc], errors="coerce").fillna(1.0)
            if p.nunique() > 1:
                fdr = multipletests(p, method="fdr_bh")[1]
                res[f"FDR_{grp}"] = fdr
            else:
                res[f"FDR_{grp}"] = np.ones(len(p))
    else:
        log.warning("⚠️ No p_* columns found — cannot compute FDRs.")

    # ensure FDR columns exist for any OR_* we will keep, even if p_* were missing
    for oc in [c for c in res.columns if c.startswith("OR_") and "_low_" not in c and "_high_" not in c]:
        grp = oc.replace("OR_", "")
        if f"FDR_{grp}" not in res.columns:
            res[f"FDR_{grp}"] = np.nan

    # log2OR
    or_cols = [c for c in res.columns if c.startswith("OR_") and "_low_" not in c and "_high_" not in c]
    for oc in or_cols:
        grp = oc.replace("OR_", "")
        x = pd.to_numeric(res[oc], errors="coerce").replace([np.inf, -np.inf], np.nan)
        res[f"log2OR_{grp}"] = np.log2(x.clip(lower=1e-30))

    keep = (["feature", "effect", "n_isolates_used", "n_patients_used", "rows",
            "mean_n", "mean_E"]
            + [c for c in res.columns if c.startswith("FDR_")]
            + [c for c in res.columns if c.startswith("log2OR_")]
            + [c for c in res.columns if c.startswith("OR_")]
            + [c for c in res.columns if c.startswith("sum_")]
            + [c for c in res.columns if c.startswith("p_")])
    res = res[[c for c in keep if c in res.columns]].copy()

    # Sort by min FDR if present
    fdr_cols = [c for c in res.columns if c.startswith("FDR_")]
    if fdr_cols:
        fdr_num = res[fdr_cols].apply(pd.to_numeric, errors="coerce")
        res["FDR_min"] = fdr_num.min(axis=1)
        res = res.sort_values(["FDR_min"] + fdr_cols, kind="mergesort")

    return res


def _to_long_from_res(res: pd.DataFrame) -> pd.DataFrame:
    # Comparisons we can actually report (must have OR_* column present)
    comps = [g for g in ("UC","CD") if f"OR_{g}" in res.columns]

    def col_or_zero(df, col):
        return (pd.to_numeric(df[col], errors="coerce").fillna(0)
                if col in df.columns else pd.Series(0, index=df.index, dtype=float))

    rows = []
    for g in [x for x in ("UC","CD") if f"OR_{x}" in res.columns]:
        ORg   = pd.to_numeric(res.get(f"OR_{g}"), errors="coerce").replace([np.inf, -np.inf], np.nan)
        FDRg  = pd.to_numeric(res.get(f"FDR_{g}"), errors="coerce")
        p_raw = pd.to_numeric(res.get(f"p_{g}"),   errors="coerce")  # <-- raw p

        rows.append(pd.DataFrame({
            "feature":     res["feature"],
            "effect":      res["effect"],
            "n_isolates":  res.get("n_isolates_used", pd.Series(np.nan, index=res.index)),
            "n_patients":  res.get("n_patients_used", pd.Series(np.nan, index=res.index)),
            "Comparison":  g,
            "log2OR":      np.log2(ORg.clip(lower=1e-30)),
            "p":           p_raw,          # <-- raw p-value
            "FDR":         FDRg,           # adjusted p
            "case_y":      pd.to_numeric(res.get(f"sum_y_{g}"), errors="coerce"),
            "case_n":      pd.to_numeric(res.get(f"sum_n_{g}"), errors="coerce"),
            "case_E":      pd.to_numeric(res.get(f"sum_E_{g}"), errors="coerce"),
            "nonIBD_y":    pd.to_numeric(res.get("sum_y_nonIBD"), errors="coerce"),
            "nonIBD_n":    pd.to_numeric(res.get("sum_n_nonIBD"), errors="coerce"),
            "nonIBD_E":    pd.to_numeric(res.get("sum_E_nonIBD"), errors="coerce"),
        }))

    out = (pd.concat(rows, ignore_index=True) if rows else
           pd.DataFrame(columns=["feature","effect","n_isolates","n_patients","Comparison",
                                 "log2OR","FDR","case_y","case_n","case_E",
                                 "nonIBD_y","nonIBD_n","nonIBD_E"]))
    out["neglog10FDR"] = -np.log10(np.clip(out["FDR"].astype(float), 1e-300, 1.0))
    return out


def main():
    """
    IBD vs nonIBD analysis using binomial Generalized Linear Models with clustered standard errors.

    For each feature (e.g., gene family, pathway), we fit a binomial GLM:
    y/n ~ group + log(E)
    where:
    - y = number of swept genes for this feature in this isolate
    - n = number of genes annotated to this feature (opportunity)
    - group = IBD subtype (UC, CD) vs nonIBD control
    - log(E) = log of total swept genes across all features (exposure)
    - Uses quasi-binomial family to handle overdispersion

    We use Jeffreys smoothing to handle separation and cluster-robust standard errors
    to account for patient-level clustering. No MAG fixed effects are included.
    """
    # ---- args & io ----
    ap = argparse.ArgumentParser(description="Compare IBD vs nonIBD per feature using binomial GLMs with clustered SEs")
    ap.add_argument("--parquet-dir", type=Path, required=True)
    ap.add_argument("--meta-file", type=Path, required=True)
    ap.add_argument("--mag-genes-file", type=Path, required=True, help="Path to all_genes.tsv file")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--filter-file", type=Path, help="Optional isolate list (iso_id OR patient_id,bin) to keep")
    ap.add_argument("--family", choices=["kegg","go","ec","gene"], default="gene")
    ap.add_argument("--min-isolates", type=int, default=10, help="min isolates with positive opportunity")
    ap.add_argument("--grouping", choices=["ibd2","dx3"], default="dx3",
                    help="ibd2 collapses UC/CD into IBD; dx3 keeps nonIBD/UC/CD")
    ap.add_argument("--background", choices=["gene", "snv_blocks"], default="gene",
                    help="Primary: 'gene' (callable genes). 'snv_blocks' collapses SNVs to LD blocks.")
    ap.add_argument("--use-E-excl", action="store_true", default=True,
                    help="Use E_excl (exposure excluding current feature) instead of E to avoid over-adjustment")
    ap.add_argument("--no-E", action="store_true", default=False,
                    help="Run sensitivity analysis without E term at all")
    ap.add_argument("--n-jobs", type=int, default=1, help="Parallel jobs for per-feature GLM fits")
    
    # Optional tuning args for SNV blocks
    ap.add_argument("--block-size-bp", type=int, default=200, help="Block size in bp for SNV blocks (default 200)")
    ap.add_argument("--max-bp-gap", type=int, default=None, help="Maximum gap in bp between SNVs in same block (default None)")
    ap.add_argument("--max-freq-diff", type=float, default=0.15, help="Maximum frequency difference between SNVs in same block (default 0.15)")
    ap.add_argument("--scaffold-file", type=Path, help="Path to combined_processed_scaffolds.txt for callable gene filtering")
    ap.add_argument("--min-breadth", type=float, default=0.80, help="Minimum breadth for callable genes (default 0.80)")
    ap.add_argument("--min-depth", type=float, default=10.0, help="Minimum coverage depth for callable genes (default 10.0)")
    ap.add_argument("--gene-len-source", choices=["snv","default"], default="snv",
                    help="Use SNV coordinates (snv) or a constant default length for all genes (default).")
    
    args = ap.parse_args()
    
    # Store args globally for parallel access
    argparse._parsed_args = args
    
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== IBD vs nonIBD Analysis ===")
    log.info("Inputs: parquet=%s meta=%s genes=%s", args.parquet_dir, args.meta_file, args.mag_genes_file)
    log.info("Family=%s  Grouping=%s  MinIsolates=%d",
             args.family, args.grouping, args.min_isolates)

    # ---- metadata ----
    log.info("[LOAD] Loading metadata from: %s", args.meta_file)
    meta_raw = read_table_auto(args.meta_file)
    meta = map_dx_to_group(meta_raw, grouping=args.grouping)
    log.info("[LOAD] Metadata loaded: %d patients, groups: %s", 
             len(meta), meta["group"].value_counts().to_dict())

    # ---- SNVs ----
    log.info("[LOAD] Loading SNVs from: %s", args.parquet_dir)
    cols_needed = ["patient_id","bin","locus_tag","freq_range","start","stop"]
    snvs = read_snvs(args.parquet_dir, cols_needed)
    if snvs.empty:
        raise SystemExit("No SNVs loaded.")
    log.info("[LOAD] SNVs loaded: %d rows, freq_range: %.3f - %.3f", 
             len(snvs), snvs["freq_range"].min(), snvs["freq_range"].max())
    
    snvs["_w"] = 1.0
    snvs = split_locus_tags(snvs, weight_col="_w")
    snvs = apply_filter(snvs, args.filter_file)
    snvs["iso_id"] = snvs["patient_id"].astype(str) + "|" + snvs["bin"].astype(str)
    log.info("[LOAD] SNVs after processing: %d rows, %d unique isolates", 
             len(snvs), snvs["iso_id"].nunique())
    
    # ---- MAG genes & feature parsing ----
    log.info("[LOAD] Loading MAG genes from: %s", args.mag_genes_file)
    mag_genes = load_mag_genes(args.mag_genes_file)
    if args.filter_file:
        mag_genes = apply_filter(mag_genes, args.filter_file)
    log.info("[LOAD] MAG genes loaded: %d rows, %d unique isolates", 
             len(mag_genes), mag_genes["iso_id"].nunique())

    # ---- Scaffold coverage & callable genes ----
    callable_genes = None
    if args.scaffold_file and args.background == "gene":
        scaffold_cov = load_scaffold_coverage(args.scaffold_file, args.meta_file)
        if args.filter_file:
            scaffold_cov = apply_filter(scaffold_cov, args.filter_file)
        
        log.info("[LOAD] Creating callable genes with min_breadth=%.2f, min_depth=%.1f", 
                 args.min_breadth, args.min_depth)
        use_snv_lens = (args.gene_len_source == "snv")
        callable_genes = make_callable_genes(
            mag_genes, scaffold_cov, snvs,
            min_breadth=args.min_breadth, 
            min_depth=args.min_depth,
            use_snv_lens=use_snv_lens
        )
        log.info("[LOAD] Callable genes created: %d gene-isolate pairs, %d callable (%.1f%%)", 
                 len(callable_genes), callable_genes["is_callable"].sum(),
                 100 * callable_genes["is_callable"].sum() / len(callable_genes) if len(callable_genes) > 0 else 0)

    # attach rep_gene/rep_dbxrefs to SNVs and parse features
    snvs_before = len(snvs)
    snvs = snvs.merge(
        mag_genes[["iso_id","locus_tag","rep_gene","rep_dbxrefs"]],
        on=["iso_id","locus_tag"], how="left"
    )
    log.info("[LOAD] SNVs after MAG gene annotation: %d -> %d rows (%.1f%% retained)", 
             snvs_before, len(snvs), 100 * len(snvs) / snvs_before if snvs_before > 0 else 0)
    
    if args.family == "gene":
        snvs["feature_list"] = snvs["rep_gene"].apply(lambda x: parse_feature_list(x, "gene"))
    else:
        snvs["feature_list"] = snvs["rep_dbxrefs"].apply(lambda x: parse_feature_list(x, args.family))
    
    snvs_before = len(snvs)
    snvs = snvs.merge(meta, on="patient_id", how="left").dropna(subset=["group"])
    log.info("[LOAD] SNVs after metadata merge and group filtering: %d -> %d rows (%.1f%% retained)", 
             snvs_before, len(snvs), 100 * len(snvs) / snvs_before if snvs_before > 0 else 0)
    
    # ---- aggregate ----
    # Choose aggregation based on background type
    if args.background == "gene":
        agg = build_agg_gene(snvs, meta, args.family, mag_genes, callable_genes=callable_genes)
    elif args.background == "snv_blocks":
        agg = build_agg_snv_blocks(
            snvs, meta, args.family, mag_genes,
            callable_genes=callable_genes,  # passes callable lengths if present
            freq_cut=FREQ_CUT,
            block_size_bp=args.block_size_bp,
            max_bp_gap=args.max_bp_gap,
            max_freq_diff=args.max_freq_diff
        )

    log.info("Aggregation rows=%d  features=%d  isolates=%d",
            len(agg), agg["feature"].nunique(), agg["iso_id"].nunique())
    
    # Diagnostic summary by group
    with log_time("[DIAG] Group summary"):
        grp = agg.groupby("group", observed=True)
        diag_summary = pd.DataFrame({
            "isolates": grp["iso_id"].nunique(),
            "sum_y":    grp["y"].sum(),
            "sum_n":    grp["n"].sum(),
            "rate":     (grp["y"].sum() / grp["n"].sum()).astype(float),
            "median_E": grp["E"].median(),
            "median_E_excl": grp["E_excl"].median() if "E_excl" in agg.columns else grp["E"].median()
        }).round(4)
        log.info("\n%s", diag_summary)

        # Top features by y share per group (quick smoke test for reversals)
        topk = (agg.assign(y_share=agg["y"]/agg["n"].clip(lower=1))
                  .sort_values("y_share", ascending=False)
                  .groupby("group", observed=True).head(5)
                  .loc[:, ["group","feature","y","n","y_share"]])
        log.info("[DIAG] Top features by y/n within each group (head):\n%s", topk.to_string(index=False))
    
    # Sanity checks
    log.info("=== SANITY CHECKS ===")
    log.info("Group categories: %s", agg["group"].cat.categories.tolist())
    log.info("Meta group counts: %s", meta["group"].value_counts().to_dict())
    log.info("Agg group counts: %s", agg["group"].value_counts().to_dict())
    
    # ---- fit, adjust, reshape ----
    res, agg_kept = _fit_all_features(agg, args.min_isolates, 
                                     use_E_excl=args.use_E_excl, no_E=args.no_E)
    log.info(f"Preview of res columns before FDR step: {res.columns.tolist()}")
    log.info(f"First 5 p-values: {res.filter(like='p_').head()}")
    res = _compute_fdrs_and_effects(res)
    long_out = _to_long_from_res(res)

    # ---- write ----
    bg_tag = getattr(args, "background", "gene")   # if you added the background arg
    group_tag = args.grouping
    fam_tag = args.family

    fname = f"ibd_vs_nonibd_{fam_tag}_{group_tag}_{bg_tag}_min{args.min_isolates}.csv"
    out_csv = out_dir / fname
    long_out.to_csv(out_csv, index=False)
    log.info("✓ Wrote %s (rows=%d)", out_csv, len(long_out))

    # Calculate phi statistics for summary
    phi_stats = {}
    if "phi" in res.columns and not res["phi"].isna().all():
        phi_vals = res["phi"].dropna()
        phi_stats = {
            "median_phi": float(phi_vals.median()),
            "mean_phi": float(phi_vals.mean()),
            "pct_overdispersed": float(100 * (phi_vals > 1.5).mean())
        }
    
    # Calculate sweep statistics by group for summary
    sweep_stats = {}
    try:
        from scipy.stats import mannwhitneyu
        
        # Calculate median sweeps per isolate by group
        group_sweeps = agg_kept.groupby(['iso_id', 'group'])['y'].sum().reset_index()
        group_medians = group_sweeps.groupby('group')['y'].median()
        
        # Calculate delta between groups
        if 'UC' in group_medians and 'nonIBD' in group_medians:
            uc_delta = group_medians['UC'] - group_medians['nonIBD']
            uc_data = group_sweeps[group_sweeps['group'] == 'UC']['y']
            nonibd_data = group_sweeps[group_sweeps['group'] == 'nonIBD']['y']
            if len(uc_data) > 0 and len(nonibd_data) > 0:
                _, uc_wilcoxon_p = mannwhitneyu(uc_data, nonibd_data, alternative='two-sided')
            else:
                uc_wilcoxon_p = np.nan
        else:
            uc_delta = np.nan
            uc_wilcoxon_p = np.nan
            
        if 'CD' in group_medians and 'nonIBD' in group_medians:
            cd_delta = group_medians['CD'] - group_medians['nonIBD']
            cd_data = group_sweeps[group_sweeps['group'] == 'CD']['y']
            nonibd_data = group_sweeps[group_sweeps['group'] == 'nonIBD']['y']
            if len(cd_data) > 0 and len(nonibd_data) > 0:
                _, cd_wilcoxon_p = mannwhitneyu(cd_data, nonibd_data, alternative='two-sided')
            else:
                cd_wilcoxon_p = np.nan
        else:
            cd_delta = np.nan
            cd_wilcoxon_p = np.nan
        
        sweep_stats = {
            "median_sweeps_per_isolate": {
                "UC": float(group_medians.get('UC', np.nan)),
                "CD": float(group_medians.get('CD', np.nan)),
                "nonIBD": float(group_medians.get('nonIBD', np.nan))
            },
            "median_delta_vs_nonIBD": {
                "UC": float(uc_delta),
                "CD": float(cd_delta)
            },
            "wilcoxon_p_vs_nonIBD": {
                "UC": float(uc_wilcoxon_p),
                "CD": float(cd_wilcoxon_p)
            }
        }
    except ImportError:
        log.warning("scipy.stats not available for sweep statistics")
        sweep_stats = {}

    with open(out_dir / "summary.json", "w") as fh:
        json.dump({
            "family": args.family,
            "engine": "binom_logit_clustered",
            "background": args.background,
            "uses_logE": True,
            "uses_quasi_binomial": True,
            "winsorization_p": 0.995,
            "parquet_dir": str(args.parquet_dir),
            "n_features_fit": int(long_out["feature"].nunique()),
            "n_isolates": int(agg_kept["iso_id"].nunique()),
            "phi_stats": phi_stats,
            "sweep_stats": sweep_stats,
            "filters": {
                "min_isolates": args.min_isolates,
                "freq_cut": FREQ_CUT,
                "min_events": MIN_EVENTS,
                "min_opp_per_group": MIN_OPP_PER_GRP,
                "each_group_in_≥2_MAGs": True
            }
        }, fh, indent=2)

if __name__ == "__main__":
    main()
