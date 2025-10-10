#!/usr/bin/env python3
"""
ibd_vs_nonibd_ko_glm.py

Fast KO/GO/EC/gene comparison between IBD and nonIBD using stratified GLMs.

Model (per feature k):
  log E[Y_ik] = α + β * IBD_i + FE(MAG_i) + log(opportunity_ik)
  where Y_ik = # sweeping variants mapped to feature k in isolate i.

Key choices:
  • Opportunity (offset): by default, the # non-sweeping variants mapped to that feature
    in the same isolate (KO-specific “exposure”). Alternatives supported:
        - ko_genes: external table of MAG×feature gene counts (best if available)
        - total_vars: total variants in isolate (coarser)
  • Stratification: MAG fixed effects via C(mag_id) in a Poisson GLM (default).
    If many strata or memory-sensitive, use --engine gee for a GEE Poisson with
    MAG clustering and robust SEs (no fixed effects).

Inputs
------
--parquet-dir        Path to parquet directory (or a single parquet file). Must contain:
                     patient_id, bin (MAG), is_sweep (bool/int), and one of:
                         kegg_terms | go_terms | ec_terms | gene
--meta-file          CSV/TSV with patient_id and diagnosis (CD/UC→IBD; nonIBD kept)
--filter-file        Optional CSV/TSV listing isolates to include (patient_id,bin or iso_id)
--opp-genes-file     Optional TSV: mag_id (bin), feature, n_genes (for --opp ko_genes)

Main options
------------
--family             one of: kegg (default), go, ec, gene
--engine             poisson_fe (default) or gee
--opp                opportunity type: non_sweep (default), ko_genes, total_vars
--min-strata         minimum # unique MAGs contributing to a feature (default 5)
--min-isolates       minimum # isolates with positive opportunity for a feature (default 20)
--sparse-cutoff      if total sweeps for a feature <= this, use binomial FE fallback (default 5)

Outputs
-------
out_dir/
  ibd_vs_nonibd_<family>.csv     # KO-level results with RR, p, FDR, CIs
  summary.json                   # run stats
  volcano_<family>.png           # log2RR vs -log10 p (optional)
"""

from __future__ import annotations
import argparse, json, logging, os, re
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import re
from scipy import stats

try:
    import pyarrow.dataset as ds
    HAVE_ARROW = True
except Exception:
    HAVE_ARROW = False

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Poisson, Binomial, Tweedie
from statsmodels.genmod.families.links import log as LogLink
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.sandwich_covariance import cov_hc0, cov_cluster_2groups, cov_cluster
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from scipy.stats import brunnermunzel
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------- logging ---------------------------------
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger("ibd_glm")

# --------------------------- feature parsing patterns ---------------------------
KO_PAT = re.compile(r"^(K\d{5})$")
GO_PAT = re.compile(r"^(\d{7})$")
EC_PAT = re.compile(r"^(\d+(?:\.\d+){0,3})$")

def _split_xrefs(dbx: object) -> List[str]:
    """Split dbxrefs string into individual tokens"""
    if pd.isna(dbx) or dbx is None:
        return []
    s = str(dbx).strip()
    if not s or s.lower() in {"nan", "none"}:
        return []
    # Split by common delimiters and clean
    tokens = []
    for part in re.split(r'[,;|\s]+', s):
        part = part.strip()
        if part and part.lower() not in {"nan", "none"}:
            tokens.append(part)
    return tokens

def feats_from_dbxrefs(dbx: object, kind: str) -> List[str]:
    """Extract features of specified kind from dbxrefs string"""
    vals = []
    for tok in _split_xrefs(dbx):
        t = tok.strip()
        # tokens can be like "KEGG:K07165", "GO:0003678", "EC:5.6.2.4", "COG:...", "RefSeq:..."
        if kind == "kegg":
            m = KO_PAT.match(t.split(":")[-1] if ":" in t else t)
            if m: 
                ko_id = "kegg:" + m.group(1).upper()
                vals.append(ko_id)
            elif t.upper().startswith("KEGG:"):
                m = KO_PAT.match(t.split(":",1)[1])
                if m: 
                    ko_id = "kegg:" + m.group(1).upper()
                    vals.append(ko_id)
        elif kind == "go":
            # Keep as lower "go:nnnnnnn"
            if t.upper().startswith("GO:"):
                m = GO_PAT.match(t.split(":",1)[1])
                if m: vals.append("go:" + m.group(1))
            else:
                m = GO_PAT.match(t)
                if m: vals.append("go:" + m.group(1))
        elif kind == "ec":
            if t.upper().startswith("EC:"):
                m = EC_PAT.match(t.split(":",1)[1])
                if m: vals.append("ec:" + m.group(1))
            else:
                m = EC_PAT.match(t)
                if m: vals.append("ec:" + m.group(1))
    return list(dict.fromkeys(vals))  # unique, keep order

def parse_feature_list(x, family: str) -> List[str]:
    """Parse feature list from various input types"""
    if x is None or (isinstance(x, float) and np.isnan(x)): 
        return []
    
    # Handle numpy arrays (common in parquet files)
    if isinstance(x, np.ndarray):
        if x.size == 0: return []
        x = list(x)
    
    # Handle lists
    if isinstance(x, list):
        if len(x) == 0: return []
        # Join list elements and split by comma
        s = ",".join(str(item) for item in x)
    else:
        s = str(x).strip()
    
    if not s: return []
    
    # For gene features, use rep_gene directly
    if family == "gene":
        if s.lower() in {"nan", "none", ""}:
            return []
        # Split by common delimiters
        toks = []
        for part in re.split(r'[,;|/\s]+', s):
            part = part.strip()
            if part and part.lower() not in {"nan", "none", ""}:
                toks.append(part)
        return toks
    
    # For other features, use dbxrefs parsing
    return feats_from_dbxrefs(x, family)


# --------------------------- helpers ---------------------------------
def read_table_auto(p: Path) -> pd.DataFrame:
    """Read CSV or TSV with minimal guessing."""
    # Check file content to determine separator
    with open(p, 'r') as f:
        first_line = f.readline().strip()
    
    # If first line contains tabs, use TSV; otherwise use CSV
    if '\t' in first_line:
        df = pd.read_csv(p, sep="\t")
        log.info(f"Loaded TSV with columns: {list(df.columns)}")
    else:
        df = pd.read_csv(p, sep=",")
        log.info(f"Loaded CSV with columns: {list(df.columns)}")
    
    return df

def clean_cols(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = (out.columns.str.strip().str.lower()
                   .str.replace(r'[^0-9a-zA-Z]+', '_', regex=True)
                   .str.replace(r'_+', '_', regex=True).str.strip('_'))
    return out

def map_dx_to_group(meta: pd.DataFrame, grouping: str = "dx3") -> pd.DataFrame:
    log.info(f"Original metadata columns: {list(meta.columns)}")
    if "patient_id" not in meta.columns:
        for c in ["Participant ID","participant_id","patient_id","Run","run","subject","sample"]:
            if c in meta.columns:
                meta = meta.rename(columns={c:"patient_id"})
                break
    if "patient_id" not in meta.columns or "diagnosis" not in meta.columns:
        raise ValueError("meta must have patient_id and diagnosis")
    
    dx = meta["diagnosis"].astype(str).str.strip().str.lower()
    if grouping == "ibd2":
        mapping = {
            "cd":"IBD","crohn":"IBD","crohn's":"IBD","crohns":"IBD",
            "uc":"IBD","ulcerative colitis":"IBD",
            "nonibd":"nonIBD","non-ibd":"nonIBD","healthy":"nonIBD","control":"nonIBD"
    }
        cats = ["nonIBD","IBD"]
    else:  # dx3
        mapping = {
            "cd":"CD","crohn":"CD","crohn's":"CD","crohns":"CD",
            "uc":"UC","ulcerative colitis":"UC",
            "nonibd":"nonIBD","non-ibd":"nonIBD","healthy":"nonIBD","control":"nonIBD"
        }
        cats = ["nonIBD","UC","CD"]
    meta["group"] = dx.map(mapping)
    meta = meta.dropna(subset=["group"])
    meta["group"] = pd.Categorical(meta["group"], categories=cats, ordered=False)
    return meta[["patient_id","group"]].drop_duplicates("patient_id", keep="first")

def build_iso_feature_presence(mag_genes: pd.DataFrame, family: str) -> pd.DataFrame:
    """
    Build iso×feature presence map from MAG genes.
    Returns DataFrame with iso_id and feature columns where each row means
    "this feature exists in this isolate's MAG".
    """
    # Feature mapping (reuse existing logic)
    if family == "gene":
        feats = mag_genes["rep_gene"].apply(lambda x: parse_feature_list(x, "gene"))
    else:
        feats = mag_genes["rep_dbxrefs"].apply(lambda x: parse_feature_list(x, family))
    
    fm = (mag_genes.assign(feature_list=feats)
                    .explode("feature_list")
                    .dropna(subset=["feature_list"]))
    fm = fm.rename(columns={"feature_list":"feature"})
    fm["iso_id"] = fm["patient_id"].astype(str) + "|" + fm["bin"].astype(str)
    present = fm[["iso_id","feature"]].drop_duplicates()
    return present

def build_patient_isolate_burden_table(snvs: pd.DataFrame, freq_col: str, present_iso_feat: pd.DataFrame) -> pd.DataFrame:
    """
    Returns one row per (patient_id, bin, group, feature) with burden and rate metrics.
    Burden = sum(freq_range) / exposure (normalized by background)
    Rate = count(variants) / exposure (pure enrichment)
    
    Only includes patient-bin-feature combinations where the feature exists in the MAG.
    Burden=0 only when feature exists in MAG but has no variants.
    
    Assumes snvs already has: patient_id, bin, group, feature_list exploded to 'feature', and freq_col in [0,1].
    """
    # Aggregate SNV signal per iso×feature
    df_features = (snvs.loc[snvs["feature_list"].str.len().astype(int) > 0,
                            ["patient_id","bin","group","feature_list", freq_col]]
                     .explode("feature_list")
                     .rename(columns={"feature_list":"feature"}))
    df_features[freq_col] = pd.to_numeric(df_features[freq_col], errors="coerce").clip(0,1)
    agg = (df_features.groupby(["patient_id","bin","group","feature"], observed=True)[freq_col]
                     .agg(freqsum="sum", n_vars="size")
                     .reset_index())

    # Attach iso_id + isolate-level exposure (CORRECTED)
    agg["iso_id"] = agg["patient_id"].astype(str) + "|" + agg["bin"].astype(str)
    tot_vars_iso = (snvs.groupby(["patient_id","bin"], observed=True)
                      .size().rename("tot_vars_iso").reset_index())
    agg = agg.merge(tot_vars_iso, on=["patient_id","bin"], how="left")
    agg["exposure"] = agg["tot_vars_iso"].astype(float).clip(lower=1.0)

    # Add zeros only where the feature is present in the MAG but no variant mapped
    iso2grp = snvs[["patient_id","bin","group"]].drop_duplicates()
    present = present_iso_feat.merge(iso2grp.assign(
                iso_id=iso2grp["patient_id"].astype(str)+"|"+iso2grp["bin"].astype(str)),
                on="iso_id", how="inner")
    present = present[["patient_id","bin","group","feature"]].drop_duplicates()

    out = present.merge(agg, on=["patient_id","bin","group","feature"], how="left")
    out[["freqsum","n_vars"]] = out[["freqsum","n_vars"]].fillna(0)
    out = out.merge(tot_vars_iso, on=["patient_id","bin"], how="left", suffixes=("","_x"))
    out["exposure"] = out["tot_vars_iso"].astype(float).clip(lower=1.0)

    out["burden"] = out["freqsum"] / out["exposure"]
    out["rate"]   = out["n_vars"] / out["exposure"]
    return out

def cliffs_delta(x, y):
    """Returns Cliff's delta in [-1,1]"""
    x = np.asarray(x); y = np.asarray(y)
    m = x.size; n = y.size
    if m == 0 or n == 0: return np.nan
    # vectorized rank trick
    # slower but simple version (fine for per-feature small Ns):
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (m * n)

def bm_two_sample(x, y):
    """Brunner–Munzel, two-sided; return stat, p"""
    # fall back to np.nan if variance is zero in both
    if (np.std(x) == 0 and np.std(y) == 0): 
        return np.nan, np.nan
    stat, p = brunnermunzel(x, y, alternative="two-sided", distribution="t")
    return float(stat), float(p)

def quantile_shift(x, y, q=0.8):
    """Difference in quantiles (y - x)"""
    return float(np.quantile(y, q) - np.quantile(x, q))

def bootstrap_ci(func, x, y, B=2000, alpha=0.05, seed=1):
    """Bootstrap confidence interval for a function of two samples"""
    rng = np.random.default_rng(seed)
    m, n = len(x), len(y)
    if m == 0 or n == 0:
        return (np.nan, np.nan)
    vals = []
    for _ in range(B):
        xb = x[rng.integers(0, m, m)]
        yb = y[rng.integers(0, n, n)]
        vals.append(func(xb, yb))
    lo, hi = np.quantile(vals, [alpha/2, 1-alpha/2])
    return float(lo), float(hi)

def run_burden_tests(burden_df: pd.DataFrame, min_patients_per_group: int = 15, 
                    min_signal_pct: float = 0.10, q: float = 0.80) -> pd.DataFrame:
    """
    Run burden and rate tests for each feature comparing groups.
    
    Parameters:
    - burden_df: DataFrame with columns [patient_id, bin, group, feature, burden, rate, exposure]
    - min_patients_per_group: minimum patients per group to test
    - min_signal_pct: minimum fraction of patients with non-zero burden in either group
    - q: quantile for tail shift calculation
    
    Returns:
    - DataFrame with test results
    """
    results = []
    
    for feat, dfk in burden_df.groupby("feature", sort=False):
        # Test burden and rate for each comparison
        for grp in [g for g in dfk["group"].unique() if g != "nonIBD"]:
            for metric in ["burden", "rate"]:
                x = dfk.loc[dfk["group"]=="nonIBD", metric].to_numpy()
                y = dfk.loc[dfk["group"]==grp, metric].to_numpy()
                
                # Hard filters
                if (len(x) < min_patients_per_group) or (len(y) < min_patients_per_group):
                    continue
                if (np.std(x)==0 and np.std(y)==0):
                    continue
                
                # Signal filter: require some non-zero values in either group
                x_nonzero = np.sum(x > 0)
                y_nonzero = np.sum(y > 0)
                if (x_nonzero / len(x) < min_signal_pct and y_nonzero / len(y) < min_signal_pct):
                    continue
                
                # Run tests
                stat, p = bm_two_sample(x, y)
                d = cliffs_delta(x, y)
                dq = quantile_shift(x, y, q=q)
                
                # Calculate prevalence and positive-only summaries (clean output)
                # Get total patient-bin combinations per group from the burden_df
                total_nonIBD = len(burden_df[(burden_df["group"] == "nonIBD") & 
                                           (burden_df["feature"] == feat)])
                total_grp = len(burden_df[(burden_df["group"] == grp) & 
                                        (burden_df["feature"] == feat)])
                
                # Positive values only (exclude zeros for meaningful medians/means)
                x_pos = x[x > 0] if len(x) > 0 else np.array([])
                y_pos = y[y > 0] if len(y) > 0 else np.array([])
                
                # Prevalence (fraction with positive values)
                prev_nonIBD = float(len(x_pos) / total_nonIBD) if total_nonIBD > 0 else 0.0
                prev_grp = float(len(y_pos) / total_grp) if total_grp > 0 else 0.0
                
                results.append({
                    "feature": feat,
                    "comparison": f"{grp}_vs_nonIBD",
                    "metric": metric,
                    "n_nonIBD_pos": int(len(x_pos)),  # patients with positive values
                    "n_"+grp+"_pos": int(len(y_pos)), # patients with positive values
                    "total_nonIBD": int(total_nonIBD),
                    "total_"+grp: int(total_grp),
                    "prev_nonIBD": prev_nonIBD,  # prevalence (fraction with positive values)
                    "prev_"+grp: prev_grp,       # prevalence (fraction with positive values)
                    "BM_stat": stat,
                    "p": p,
                    "cliffs_delta": d,
                    f"delta_Q{int(100*q)}": dq,
                    # Positive-only summaries (meaningful medians/means)
                    "median_nonIBD_pos": float(np.median(x_pos)) if len(x_pos) > 0 else np.nan,
                    "median_"+grp+"_pos": float(np.median(y_pos)) if len(y_pos) > 0 else np.nan,
                    "mean_nonIBD_pos": float(np.mean(x_pos)) if len(x_pos) > 0 else np.nan,
                    "mean_"+grp+"_pos": float(np.mean(y_pos)) if len(y_pos) > 0 else np.nan,
                    # All-inclusive summaries (including zeros) - for reference
                    "median_nonIBD_all": float(np.median(x)) if len(x) > 0 else np.nan,
                    "median_"+grp+"_all": float(np.median(y)) if len(y) > 0 else np.nan,
                    "mean_nonIBD_all": float(np.mean(x)) if len(x) > 0 else np.nan,
                    "mean_"+grp+"_all": float(np.mean(y)) if len(y) > 0 else np.nan,
                })
    
    return pd.DataFrame(results)

def build_variant_table(snvs: pd.DataFrame, freq_col: str) -> pd.DataFrame:
    """
    One row per (iso_id, mag_id, group, feature, locus_tag) with a single
    representative freq (mean if multiple variants map to the same locus in that isolate).
    Ensures freq in [0,1].
    """
    req = {"patient_id","bin","locus_tag","feature_list","mag_id","group",freq_col}
    missing = req - set(snvs.columns)
    if missing:
        raise ValueError(f"build_variant_table missing columns: {missing}")

    df = snvs.loc[snvs["feature_list"].str.len().astype(int) > 0,
                  ["patient_id","bin","locus_tag","feature_list","mag_id","group",freq_col]].copy()
    df["iso_id"] = df["patient_id"].astype(str) + "|" + df["bin"].astype(str)

    df = df.explode("feature_list").rename(columns={"feature_list":"feature"})
    df[freq_col] = pd.to_numeric(df[freq_col], errors="coerce").clip(0,1)
    df = df.dropna(subset=[freq_col,"feature","locus_tag","iso_id","mag_id","group"])

    # dedup within (iso, feature, locus): average
    out = (df.groupby(["iso_id","mag_id","group","feature","locus_tag"], observed=True)[freq_col]
             .mean().rename("freq").reset_index())
    return out

def fit_feature_logit_lm(dfv: pd.DataFrame) -> Optional[dict]:
    # mean shift on logit(freq) controlling for MAG FE
    if dfv["group"].nunique() < 2 or dfv["mag_id"].nunique() < 2: 
        return None
    d = dfv.copy()
    eps = 1e-6
    d["y"] = np.log((d["freq"] + eps) / (1 - d["freq"] + eps))
    # Ensure nonIBD is baseline if present
    if isinstance(d["group"].dtype, pd.CategoricalDtype):
        cats = list(d["group"].dtype.categories)
        if "nonIBD" in cats and cats[0] != "nonIBD":
            cats.remove("nonIBD"); cats = ["nonIBD"] + cats
            d["group"] = d["group"].cat.set_categories(cats, ordered=False)
    try:
        m = smf.ols("y ~ C(group) + C(mag_id)", data=d).fit()
        
        # Two-way clustering by iso_id and mag_id
        try:
            from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
            cov = cov_cluster_2groups(m, d["iso_id"], d["mag_id"])
        except:
            from statsmodels.stats.sandwich_covariance import cov_cluster
            cov = cov_cluster(m, d["iso_id"])
        
        se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
        out = {"effect":"logit_LM","n_iso":int(d["iso_id"].nunique()),
               "n_mag":int(d["mag_id"].nunique()),"n_vars":len(d)}
        base = d["group"].cat.categories[0] if hasattr(d["group"].dtype,"categories") else d["group"].unique()[0]
        for lev in [c for c in d["group"].unique() if c != base]:
            pname = f"C(group)[T.{lev}]"
            if pname in names:
                beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                z = beta / max(se, 1e-12)
                p = float(2 * (1 - stats.norm.cdf(abs(z))))
                out.update({f"beta_{lev}":beta, f"se_{lev}":se, f"p_{lev}":p})
        return out
    except Exception:
        return None

def fit_feature_qreg(dfv: pd.DataFrame, tau: float) -> Optional[dict]:
    # upper-tail shift at quantile tau, controlling for MAG FE
    if dfv["group"].nunique() < 2 or dfv["mag_id"].nunique() < 2:
        return None
    try:
        mod = smf.quantreg("freq ~ C(group) + C(mag_id)", data=dfv)
        res = mod.fit(q=tau)
        out = {"effect":f"qreg{int(100*tau)}","n_iso":int(dfv["iso_id"].nunique()),
               "n_mag":int(dfv["mag_id"].nunique()),"n_vars":len(dfv)}
        for pname in res.params.index:
            if pname.startswith("C(group)[T."):
                lev = pname.split("[T.",1)[1][:-1]
                out.update({f"beta_{lev}": float(res.params[pname]),
                            f"p_{lev}": float(res.pvalues[pname])})
        return out
    except Exception:
        return None

# Removed duplicate parse_feature_list function - using the robust version that calls feats_from_dbxrefs

def read_snvs(parquet_dir: Path, cols_needed: List[str]) -> pd.DataFrame:
    if parquet_dir.is_file() and parquet_dir.suffix == ".parquet":
        df = pd.read_parquet(parquet_dir, columns=[c for c in cols_needed if c != "*"])
        return df
    if HAVE_ARROW:
        dataset = ds.dataset(str(parquet_dir), format="parquet")
        # select columns if present
        schema_cols = set([n for n in dataset.schema.names])
        cols = [c for c in cols_needed if c in schema_cols or c == "*"]
        tbl = dataset.to_table(columns=[c for c in cols if c != "*"])
        df = tbl.to_pandas(split_blocks=True, self_destruct=True)
        return df
    # fallback: glob
    dfs = []
    for p in sorted(parquet_dir.glob("*.parquet")):
        dfs.append(pd.read_parquet(p, columns=[c for c in cols_needed if c != "*"]))
    if not dfs:
        raise ValueError(f"No parquet files in {parquet_dir}")
    return pd.concat(dfs, ignore_index=True)

def split_locus_tags(
    df: pd.DataFrame,
    col: str = "locus_tag",
    sep: str = "/",
    drop_dupes: bool = False,
    split_weight_col: str | None = None,  # if you want to divide a weight across splits
) -> pd.DataFrame:
    """
    Vectorized split of delimited locus tags into separate rows.
    This is critical for proper annotation mapping.

    Parameters
    ----------
    df : DataFrame
    col : column to split (default 'locus_tag')
    sep : delimiter (default '/')
    drop_dupes : if True, drop duplicate rows after splitting
    split_weight_col : if provided, divide this column by the number of tags
                      in the original row (keeps index to do it correctly)

    Returns
    -------
    DataFrame with one locus_tag per row, trimmed, non-empty.
    """

    if col not in df.columns:
        raise KeyError(f"Column '{col}' not in DataFrame")

    n0 = len(df)
    log.info("Splitting %s: input rows=%d", col, n0)

    # Ensure string dtype and split to lists
    s = df[col].astype("string")
    has_sep = s.str.contains(sep, na=False)
    log.info("Rows containing '%s': %d", sep, int(has_sep.sum()))

    # Precompute split counts per original row (align by index)
    split_lists = s.str.split(sep)
    split_counts = split_lists.str.len()  # Simple count of split elements

    # Join the split lists, then explode without resetting index
    out = df.join(split_lists.rename("_split_tags"))
    out = out.explode("_split_tags")  # keeps original index by default
    out[col] = out["_split_tags"].astype("string").str.strip()
    out = out.drop(columns=["_split_tags"])

    # Keep only non-empty/non-null tokens
    mask = out[col].notna() & (out[col] != "") & (~out[col].str.lower().isin({"nan", "none"}))
    out = out[mask]

    # Optionally divide a weight column across the number of tokens in the original row
    if split_weight_col is not None:
        if split_weight_col not in out.columns:
            raise KeyError(f"split_weight_col '{split_weight_col}' not in DataFrame")
        # Align counts by original index; broadcast via .loc on the current index
        out["_split_n"] = split_counts.loc[out.index].clip(lower=1)
        out[split_weight_col] = out[split_weight_col] / out["_split_n"]
        out = out.drop(columns=["_split_n"])

    # Optional deduplication (only if you truly want to collapse exact row duplicates)
    if drop_dupes:
        out = out.drop_duplicates()

    n1 = len(out)
    log.info("After split: output rows=%d (Δ=%+d)", n1, n1 - n0)
    return out

def load_mag_genes(mag_path: Path) -> pd.DataFrame:
    """Load MAG genes file and prepare for feature mapping"""
    log.info(f"Loading MAG genes from: {mag_path}")
    cols = ["patient_id","bin","scaffold","locus_tag","cluster_rep","rep_gene","rep_dbxrefs"]
    mg = pd.read_csv(mag_path, sep="\t")
    log.info(f"Loaded MAG genes. Shape: {mg.shape}")
    log.info(f"MAG gene columns: {list(mg.columns)}")
    
    missing = [c for c in cols if c not in mg.columns]
    if missing:
        raise ValueError(f"MAG gene file missing columns: {missing}")
    
    for c in ["patient_id","bin","locus_tag","cluster_rep","rep_gene","rep_dbxrefs"]:
        mg[c] = mg[c].astype(str).str.strip()
    mg["iso_id"] = mg["patient_id"] + "|" + mg["bin"]
    
    log.info(f"Unique iso_ids in MAG genes: {mg['iso_id'].nunique()}")    
    return mg[["patient_id","bin","iso_id","locus_tag","cluster_rep","rep_gene","rep_dbxrefs"]]

def apply_filter(snvs: pd.DataFrame, filter_file: Optional[Path]) -> pd.DataFrame:
    if not filter_file: return snvs
    filt = read_table_auto(filter_file)
    # Don't use clean_cols on filter - it corrupts column names
    # filt = clean_cols(filt)
    if "iso_id" in filt.columns:
        filt = filt[["iso_id"]].dropna().drop_duplicates()
        snvs["iso_id"] = snvs["patient_id"].astype(str) + "|" + snvs["bin"].astype(str)
        return snvs.merge(filt, on="iso_id", how="inner").drop(columns=["iso_id"])
    # expect patient_id, bin
    if not {"patient_id","bin"}.issubset(filt.columns):
        # try to split combined
        for col in filt.columns:
            if "patient" in col and "bin" in col:
                parts = filt[col].astype(str).str.split(r"[|,\t]", n=1, expand=True)
                filt = pd.DataFrame({"patient_id": parts[0], "bin": parts[1]})
                break
    if not {"patient_id","bin"}.issubset(filt.columns):
        raise ValueError("filter_file must have iso_id or patient_id,bin")
    filt["patient_id"] = filt["patient_id"].astype(str).str.strip()
    filt["bin"] = filt["bin"].astype(str).str.strip()
    snvs["patient_id"] = snvs["patient_id"].astype(str)
    snvs["bin"] = snvs["bin"].astype(str)
    before = len(snvs)
    out = snvs.merge(filt[["patient_id","bin"]].drop_duplicates(), on=["patient_id","bin"], how="inner")
    log.info("Applied isolate filter: %d → %d rows (%.1f%% kept)", before, len(out), 100*len(out)/max(1,before))
    return out

def build_agg(snvs: pd.DataFrame, meta: pd.DataFrame, family: str,
              opp_type: str, opp_gene_df: Optional[pd.DataFrame], 
              y_metric: str = "hits", freq_col: str = "freq_range") -> pd.DataFrame:
    """Return per-isolate×feature table with sweep_hits, opportunity, and strata."""
    # Don't use clean_cols on SNVs - it corrupts column names
    # snvs = clean_cols(snvs)
    # columns harmonization
    if "is_sweep" not in snvs.columns:
        raise ValueError("SNVs must have is_sweep")
    for need in ("patient_id","bin"):
        if need not in snvs.columns:
            raise ValueError(f"SNVs missing required column {need}")

    # add iso_id, mag_id (group already added in main)
    snvs["iso_id"] = snvs["patient_id"].astype(str) + "|" + snvs["bin"].astype(str)
    snvs["mag_id"] = snvs["bin"].astype(str)

    # Process features from MAG genes
    log.info(f"Processing {family} features from MAG genes...")
    log.info(f"SNVs shape before feature parsing: {snvs.shape}")
    
    # Filter SNVs that have features
    has_feat = snvs["feature_list"].str.len().astype(int) > 0
    log.info(f"SNVs with features: {has_feat.sum()} out of {len(has_feat)}")
    
    # Select columns based on y_metric
    feat_cols = ["iso_id","mag_id","group","is_sweep","feature_list"]
    if y_metric == "freqsum":
        if freq_col not in snvs.columns:
            raise ValueError(f"--y-metric=freqsum requires column '{freq_col}' in SNVs")
        feat_cols.append(freq_col)
    
    snvs_feat = snvs.loc[has_feat, feat_cols].copy()
    log.info(f"SNVs with features shape: {snvs_feat.shape}")
    
    exploded = (snvs_feat.explode("feature_list")
                          .rename(columns={"feature_list":"feature"}))
    log.info(f"Exploded shape: {exploded.shape}")
    
    # Handle frequency column if needed
    if y_metric == "freqsum":
        exploded[freq_col] = pd.to_numeric(exploded[freq_col], errors="coerce").fillna(0.0)
    
    # Aggregate counts per isolate×feature for sweeps and non-sweeps
    exploded["is_sweep"] = exploded["is_sweep"].astype(bool)
    grp = exploded.groupby(["iso_id","mag_id","group","feature"], observed=True)

    if y_metric == "freqsum":
        agg = grp.agg(
            sweep_hits=("is_sweep","sum"),         # keep for reference
            n=("is_sweep","size"),
            freq_sum=(freq_col,"sum"),
            freq_max=(freq_col,"max"),
            freq_mean=(freq_col,"mean"),
        ).reset_index()
    else:
        agg = grp.agg(sweep_hits=("is_sweep","sum"), n=("is_sweep","size")).reset_index()
    
    # derive non-sweep counts (opportunity proxy)
    agg["non_sweep_hits"] = agg["n"] - agg["sweep_hits"]

    # attach KO-gene opportunity if requested
    if opp_type == "ko_genes":
        if opp_gene_df is None:
            raise ValueError("--opp ko_genes requires --opp-genes-file")
        # expect columns: mag_id (bin), feature, n_genes
        # harmonize names
        col_map = {}
        if "bin" in opp_gene_df.columns: col_map["bin"] = "mag_id"
        if "ko" in opp_gene_df.columns:  col_map["ko"] = "feature"
        if "count" in opp_gene_df.columns: col_map["count"] = "n_genes"
        opp_gene_df = opp_gene_df.rename(columns=col_map)
        need = {"mag_id","feature","n_genes"}
        if not need.issubset(opp_gene_df.columns):
            raise ValueError("opp-genes-file must have mag_id, feature, n_genes")
        opp_gene_df["feature"] = opp_gene_df["feature"].astype(str)
        agg = agg.merge(opp_gene_df[["mag_id","feature","n_genes"]],
                        on=["mag_id","feature"], how="left")
        agg["opportunity"] = agg["n_genes"].fillna(0).astype(float)
    elif opp_type == "non_sweep":
        agg["opportunity"] = agg["non_sweep_hits"].astype(float) if y_metric == "hits" else agg["n"].astype(float)
    elif opp_type == "total_vars":
        # compute total variants per isolate as exposure (crude)
        tot = exploded.groupby("iso_id", observed=True).size().rename("tot_vars").reset_index()
        agg = agg.merge(tot, on="iso_id", how="left")
        agg["opportunity"] = agg["tot_vars"].astype(float)
    else:
        raise ValueError("--opp must be one of non_sweep, ko_genes, total_vars")

    # keep rows with positive opportunity or positive signal
    agg["opportunity"] = agg["opportunity"].fillna(0.0)
    if y_metric == "freqsum":
        agg = agg[(agg["opportunity"] > 0) | (agg["freq_sum"] > 0)].copy()
    else:
        agg = agg[(agg["opportunity"] > 0) | (agg["sweep_hits"] > 0)].copy()

    return agg

def fit_one_feature(dfk: pd.DataFrame, engine: str,
                    sparse_cutoff: int = 5, y_metric: str = "hits", 
                    robust_se: str = "hc0", use_gee: bool = False,
                    tweedie_p: float = 1.5) -> Optional[dict]:
    """
    Fit per-feature model and return dict with results.
    dfk columns: iso_id, mag_id, group (IBD/nonIBD), sweep_hits, opportunity
    """
    log.debug(f"Fitting feature with {len(dfk)} observations, {dfk['iso_id'].nunique()} isolates, {dfk['mag_id'].nunique()} MAGs")
    
    # groups
    glevels = dfk["group"].dropna().unique().tolist()
    if len(glevels) < 2:
        log.debug(f"Insufficient groups: {glevels}")
        return None
    
    # ensure factor coding with nonIBD as baseline
    dfk = dfk.copy()
    if not pd.api.types.is_categorical_dtype(dfk["group"]):
        cats = ["nonIBD"] + [c for c in pd.unique(dfk["group"]) if c != "nonIBD"]
        dfk["group"] = pd.Categorical(dfk["group"], categories=cats, ordered=False)
    else:
        cats = list(dfk["group"].cat.categories)
        if "nonIBD" in cats and cats[0] != "nonIBD":
            cats.remove("nonIBD"); cats = ["nonIBD"] + cats
            dfk["group"] = dfk["group"].cat.set_categories(cats, ordered=False)

    total_sweeps = int(dfk["sweep_hits"].sum())
    log.debug(f"Total sweeps: {total_sweeps}, y_metric: {y_metric}, robust_se: {robust_se}")

    # Avoid zero offset: clip very small
    eps = 1e-8
    dfk["offset_log"] = np.log(np.clip(dfk["opportunity"].astype(float), eps, None))

    # Tweedie GLM for freqsum metric
    if y_metric == "freqsum":
        log.debug("Fitting Tweedie GLM for freqsum metric")
        # Tweedie (1<p<2) handles many zeros + continuous positive
        var_power = tweedie_p
        fam = Tweedie(var_power=var_power, link=LogLink())
        
        # Make nonIBD baseline if present
        if isinstance(dfk["group"].dtype, pd.CategoricalDtype):
            cats = list(dfk["group"].dtype.categories)
            if "nonIBD" in cats and cats[0] != "nonIBD":
                cats.remove("nonIBD"); cats = ["nonIBD"] + cats
                dfk["group"] = dfk["group"].cat.set_categories(cats, ordered=False)

        # 1) Try FE after dropping singleton MAGs
        df_fe = dfk[dfk.groupby("mag_id")["iso_id"].transform("nunique") >= 2].copy()
        try:
            m = smf.glm("freq_sum ~ C(group) + C(mag_id)", data=df_fe, family=fam,
                        offset=df_fe["offset_log"]).fit(disp=0)
            # If saturated, force fallback
            if not np.isfinite(m.df_resid) or m.df_resid <= 0:
                raise ValueError("FE model saturated (df_resid<=0)")
            # robust cov (fallback to model cov if needed)
            try:
                cov = cov_hc0(m) if robust_se == "hc0" else (
                      cov_cluster(m, df_fe["mag_id"].astype(str)) if robust_se == "cluster" else m.cov_params())
            except Exception:
                cov = m.cov_params()
            se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
            out = {"effect":"tweedie_FE","n_iso":int(df_fe["iso_id"].nunique()),
                   "n_mag":int(df_fe["mag_id"].nunique()), "total_freqsum":float(dfk["freq_sum"].sum())}
            for pname in names:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                    z = beta / max(se, 1e-12)
                    p  = float(2*(1 - stats.norm.cdf(abs(z))))
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out

        except Exception:
            # 2) Fallback: no FE, cluster by MAG
            log.debug("Tweedie FE failed, using cluster fallback")
            m = smf.glm("freq_sum ~ C(group)", data=dfk, family=fam,
                        offset=dfk["offset_log"]).fit(disp=0)
            try:
                cov = cov_cluster(m, dfk["mag_id"].astype(str))
            except Exception:
                cov = m.cov_params()
            se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
            out = {"effect":"tweedie_clusterMAG","n_iso":int(dfk["iso_id"].nunique()),
                   "n_mag":int(dfk["mag_id"].nunique()), "total_freqsum":float(dfk["freq_sum"].sum())}
            for pname in names:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                    z = beta / max(se, 1e-12)
                    p  = float(2*(1 - stats.norm.cdf(abs(z))))
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out

    # Sparse fallback: Binomial FE on presence (no offset)
    if total_sweeps <= sparse_cutoff:
        log.debug(f"Using binomial FE for sparse data (total_sweeps={total_sweeps} <= {sparse_cutoff})")
        dfk["present"] = (dfk["sweep_hits"] > 0).astype(int)
        # 1) Try FE after dropping singleton MAGs
        df_fe = dfk[dfk.groupby("mag_id")["iso_id"].transform("nunique") >= 2].copy()
        try:
            m = smf.glm("present ~ C(group) + C(mag_id)",
                        data=df_fe, family=Binomial()).fit(disp=0)
            # If saturated, force fallback
            if not np.isfinite(m.df_resid) or m.df_resid <= 0:
                raise ValueError("FE model saturated (df_resid<=0)")
            
            # choose covariance
            try:
                cov = cov_hc0(m) if robust_se == "hc0" else (
                      cov_cluster(m, df_fe["mag_id"].astype(str)) if robust_se == "cluster" else m.cov_params())
            except Exception:
                cov = m.cov_params()
            se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
            out = {"effect": "binomial_FE", "n_iso": int(df_fe["iso_id"].nunique()),
                   "n_mag": int(df_fe["mag_id"].nunique()), "total_sweeps": total_sweeps}

            # pairwise vs baseline for every non-baseline level
            for pname in names:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                    z = beta / max(se, 1e-12)
                    p  = float(2*(1 - stats.norm.cdf(abs(z))))
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out

        except (Exception, PerfectSeparationError):
            # 2) Fallback: no FE, cluster by MAG
            log.debug("Binomial FE failed, using cluster fallback")
            try:
                m = smf.glm("present ~ C(group)",
                            data=dfk, family=Binomial()).fit(disp=0)
                try:
                    cov = cov_cluster(m, dfk["mag_id"].astype(str))
                except Exception:
                    cov = m.cov_params()
                se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
                out = {"effect": "binomial_clusterMAG", "n_iso": int(dfk["iso_id"].nunique()),
                    "n_mag": int(dfk["mag_id"].nunique()), "total_sweeps": total_sweeps}

                # pairwise vs baseline for every non-baseline level
                for pname in names:
                    if pname.startswith("C(group)[T."):
                        lev = pname.split("[T.",1)[1][:-1]
                        beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                        z = beta / max(se, 1e-12)
                        p  = float(2*(1 - stats.norm.cdf(abs(z))))
                        RR = float(np.exp(beta))
                        out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                    f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                    f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
                return out
            except (Exception, PerfectSeparationError):
                # 3) Final fallback: check for perfect separation and return None if so
                log.debug("Binomial cluster failed, checking for perfect separation")
                # Check if any group has all zeros or all ones
                group_stats = dfk.groupby("group")["present"].agg(["sum", "count"]).reset_index()
                group_stats["all_zeros"] = group_stats["sum"] == 0
                group_stats["all_ones"] = group_stats["sum"] == group_stats["count"]
                if group_stats["all_zeros"].any() or group_stats["all_ones"].any():
                    log.debug("Perfect separation detected, skipping feature")
                    return None
                else:
                    log.debug("Unknown error in binomial model")
            return None

    # Poisson FE (default)
    if engine == "poisson_fe":
        log.debug("Fitting Poisson GLM with fixed effects")
        # 1) Try FE after dropping singleton MAGs
        df_fe = dfk[dfk.groupby("mag_id")["iso_id"].transform("nunique") >= 2].copy()
        try:
            m = smf.glm("sweep_hits ~ C(group) + C(mag_id)",
                        data=df_fe,
                        family=Poisson(), offset=df_fe["offset_log"]).fit(disp=0)
            # If saturated, force fallback
            if not np.isfinite(m.df_resid) or m.df_resid <= 0:
                raise ValueError("FE model saturated (df_resid<=0)")
            
            # choose covariance
            try:
                cov = cov_hc0(m) if robust_se == "hc0" else (
                      cov_cluster(m, df_fe["mag_id"].astype(str)) if robust_se == "cluster" else m.cov_params())
            except Exception:
                cov = m.cov_params()
            se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
            out = {"effect": "poisson_FE", "n_iso": int(df_fe["iso_id"].nunique()),
                   "n_mag": int(df_fe["mag_id"].nunique()), "total_sweeps": total_sweeps}

            # pairwise vs baseline for every non-baseline level
            for pname in names:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                    z = beta / max(se, 1e-12)
                    p  = float(2*(1 - stats.norm.cdf(abs(z))))
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out

        except Exception:
            # 2) Fallback: no FE, cluster by MAG
            log.debug("Poisson FE failed, using cluster fallback")
            m = smf.glm("sweep_hits ~ C(group)",
                        data=dfk,
                        family=Poisson(), offset=dfk["offset_log"]).fit(disp=0)
            try:
                cov = cov_cluster(m, dfk["mag_id"].astype(str))
            except Exception:
                cov = m.cov_params()
            se_vec = np.sqrt(np.diag(cov)); names = list(m.params.index)
            out = {"effect": "poisson_clusterMAG", "n_iso": int(dfk["iso_id"].nunique()),
                    "n_mag": int(dfk["mag_id"].nunique()), "total_sweeps": total_sweeps}

            # pairwise vs baseline for every non-baseline level
            for pname in names:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(m.params[pname]); se = float(se_vec[names.index(pname)])
                    z = beta / max(se, 1e-12)
                    p  = float(2*(1 - stats.norm.cdf(abs(z))))
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":p,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out

    # GEE Poisson with MAG clustering (robust SEs)
    if engine == "gee":
        try:
            gee = smf.gee("sweep_hits ~ C(group)",
                          groups="mag_id",
                          data=dfk,
                          family=Poisson(),
                          offset=dfk["offset_log"]).fit()
            
            out = {"effect": "gee_poisson", "n_iso": int(dfk["iso_id"].nunique()),
                    "n_mag": int(dfk["mag_id"].nunique()), "total_sweeps": total_sweeps}

            # pairwise vs baseline for every non-baseline level
            for pname in gee.params.index:
                if pname.startswith("C(group)[T."):
                    lev = pname.split("[T.",1)[1][:-1]
                    beta = float(gee.params[pname])
                    se = float(gee.bse[pname])
                    pval = float(gee.pvalues[pname])
                    RR = float(np.exp(beta))
                    out.update({f"RR_{lev}":RR, f"p_{lev}":pval,
                                f"ci_low_{lev}":float(np.exp(beta - 1.96*se)),
                                f"ci_high_{lev}":float(np.exp(beta + 1.96*se))})
            return out
        except Exception as e:
            return None

    return None

def volcano(df: pd.DataFrame, out_png: Path, title: str):
    if df.empty:
        return
    d = df.replace([np.inf,-np.inf], np.nan).dropna(subset=["log2RR","p"])
    if d.empty:
        return
    plt.figure(figsize=(7.5,5.5))
    plt.scatter(d["log2RR"], -np.log10(d["p"]), s=10, alpha=0.7)
    sig = d["FDR"] < 0.05
    if sig.any():
        plt.scatter(d.loc[sig,"log2RR"], -np.log10(d.loc[sig,"p"]), s=12, alpha=0.9, color="red")
    plt.axhline(-np.log10(0.05), ls="--", lw=0.8, color="gray")
    plt.xlabel("log2 Rate Ratio (IBD / nonIBD)")
    plt.ylabel("-log10 p-value")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# --------------------------- main ---------------------------------
def main():
    ap = argparse.ArgumentParser(description="Compare IBD vs nonIBD per feature using stratified GLMs")
    ap.add_argument("--parquet-dir", type=Path, required=True)
    ap.add_argument("--meta-file", type=Path, required=True)
    ap.add_argument("--mag-genes-file", type=Path, required=True, help="Path to all_genes.tsv file")
    ap.add_argument("--out-dir", type=Path, required=True)

    ap.add_argument("--filter-file", type=Path, help="Optional isolate list (iso_id OR patient_id,bin) to keep")
    ap.add_argument("--family", choices=["kegg","go","ec","gene"], default="kegg")
    ap.add_argument("--engine", choices=["poisson_fe","gee"], default="poisson_fe",
                    help="Poisson with MAG fixed effects (default) or GEE Poisson clustering by MAG")
    ap.add_argument("--opp", choices=["non_sweep","ko_genes","total_vars"], default="non_sweep",
                    help="Opportunity definition for offset")
    ap.add_argument("--opp-genes-file", type=Path, help="TSV with mag_id/bin, feature, n_genes (for --opp ko_genes)")

    ap.add_argument("--min-strata", type=int, default=5, help="min unique MAGs contributing to a feature")
    ap.add_argument("--min-isolates", type=int, default=20, help="min isolates with positive opportunity")
    ap.add_argument("--sparse-cutoff", type=int, default=5, help="<= uses binomial FE on presence")

    ap.add_argument("--volcano", action="store_true", help="emit a volcano plot")
    ap.add_argument("--y-metric", choices=["hits","freqsum"], default="hits",
                    help="Outcome: sweep hit counts (hits) or sum of freq ranges (freqsum)")
    ap.add_argument("--freq-col", default="freq_range",
                    help="Per-variant column in parquet with frequency range in [0,1]")
    ap.add_argument("--robust-se", choices=["none","hc0","cluster"], default="hc0",
                    help="Robust standard errors: none, hc0 (heteroskedasticity-consistent), or cluster by MAG")
    ap.add_argument("--gee-engine", action="store_true",
                    help="Use GEE instead of GLM for very large MAG counts")
    
    # New analysis options
    ap.add_argument("--analysis", choices=["glm","qc"], default="glm",
                    help="glm = Main GLM analysis (default); qc = QC mode with burden distribution tests")
    ap.add_argument("--grouping", choices=["ibd2","dx3"], default="dx3",
                    help="ibd2 collapses UC/CD into IBD; dx3 keeps nonIBD/UC/CD")
    ap.add_argument("--tweedie-p", type=float, default=1.5,
                    help="Tweedie variance power for --y-metric=freqsum (1<p<2)")
    ap.add_argument("--dist-test", choices=["logit_lm","qreg80","both"], default="both",
                    help="Which distribution test(s) to run when --analysis=dist")
    ap.add_argument("--tau", type=float, default=0.80,
                    help="Quantile for --dist-test=qreg80 (e.g., 0.8 for 80th percentile)")
    ap.add_argument("--min-patients-per-group", type=int, default=15,
                    help="Minimum patients per group for dist_simple analysis")
    ap.add_argument("--q", type=float, default=0.80,
                    help="Quantile for tail shift (ΔQ) in dist_simple analysis")
    ap.add_argument("--min-signal-pct", type=float, default=0.10,
                    help="Minimum fraction of patients with non-zero burden in either group")
    
    args = ap.parse_args()
    
    # Handle GEE engine flag
    if args.gee_engine:
        args.engine = "gee"

    out_dir = args.out_dir; out_dir.mkdir(parents=True, exist_ok=True)
    
    log.info("=== IBD vs nonIBD Analysis Started ===")
    log.info(f"Input parquet directory: {args.parquet_dir}")
    log.info(f"Metadata file: {args.meta_file}")
    log.info(f"MAG genes file: {args.mag_genes_file}")
    log.info(f"Output directory: {args.out_dir}")
    log.info(f"Feature family: {args.family}")
    log.info(f"Y-metric: {args.y_metric}")
    log.info(f"Robust SE: {args.robust_se}")
    log.info(f"Engine: {args.engine}")
    log.info(f"Opportunity: {args.opp}")
    log.info(f"Min strata: {args.min_strata}, Min isolates: {args.min_isolates}")
    
    # Warn about potentially biased offset choice
    if args.y_metric == "freqsum" and args.opp == "non_sweep":
        log.warning("Using --y-metric=freqsum with --opp=non_sweep may be biased. Consider --opp=total_vars or --opp=ko_genes.")

    # Load meta and map groups
    log.info("Loading and processing metadata...")
    meta = read_table_auto(args.meta_file)
    meta = map_dx_to_group(meta, grouping=args.grouping)
    log.info("Meta: %d patients; groups: %s", meta["patient_id"].nunique(),
             meta["group"].value_counts().to_dict())

    # Read SNVs minimal columns (no feature columns needed - we'll get them from all_genes.tsv)
    log.info("Loading SNV data...")
    cols_needed = ["patient_id","bin","is_sweep","locus_tag"]
    if args.y_metric == "freqsum" or args.analysis in ["dist", "qc"]:
        cols_needed.append(args.freq_col)
        log.info(f"Adding frequency column: {args.freq_col}")
    log.info(f"Required SNV columns: {cols_needed}")
    
    snvs = read_snvs(args.parquet_dir, cols_needed)
    if snvs.empty: raise SystemExit("No SNVs loaded.")
    log.info("SNVs loaded: %d rows", len(snvs))
    log.info(f"SNV columns: {list(snvs.columns)}")

    # CRITICAL: Split locus tags containing '/' into separate rows
    log.info("Splitting locus tags...")
    snvs = split_locus_tags(snvs)
    
    # Optional filter
    log.info(f"Applying isolate filter: {args.filter_file}")
    snvs = apply_filter(snvs, args.filter_file)
    
    # Add iso_id to SNVs before merging
    snvs["iso_id"] = snvs["patient_id"].astype(str) + "|" + snvs["bin"].astype(str)
    
    # Load MAG genes for feature mapping
    log.info("Loading MAG genes for feature mapping...")
    mag_genes = load_mag_genes(args.mag_genes_file)
    log.info(f"MAG genes loaded: {len(mag_genes)} rows")
    log.info(f"MAG genes columns: {list(mag_genes.columns)}")
    
    # Apply same filter to MAG genes
    if args.filter_file:
        log.info("Applying same filter to MAG genes...")
        mag_genes = apply_filter(mag_genes, args.filter_file)
        log.info(f"MAG genes after filtering: {len(mag_genes)} rows")
    
    # Merge SNVs with MAG genes to get feature information
    log.info("Merging SNVs with MAG genes for feature mapping...")
    log.info(f"SNVs before merge: {len(snvs)} rows")
    log.info(f"MAG genes for merge: {len(mag_genes)} rows")
    
    snvs_with_features = snvs.merge(
        mag_genes[["iso_id", "locus_tag", "rep_gene", "rep_dbxrefs"]], 
        on=["iso_id", "locus_tag"], 
        how="left"
    )
    log.info(f"SNVs after merge: {len(snvs_with_features)} rows")
    log.info(f"Merged columns: {list(snvs_with_features.columns)}")
    
    # Parse features based on family type
    log.info(f"Parsing {args.family} features from MAG genes...")
    if args.family == "gene":
        snvs_with_features["feature_list"] = snvs_with_features["rep_gene"].apply(
            lambda x: parse_feature_list(x, args.family)
        )
    else:
        snvs_with_features["feature_list"] = snvs_with_features["rep_dbxrefs"].apply(
            lambda x: parse_feature_list(x, args.family)
        )
    
    snvs = snvs_with_features
    log.info(f"Final SNVs dataset: {len(snvs)} rows")

    # Add group and mag_id to SNVs for both analysis types
    snvs = snvs.merge(meta, on="patient_id", how="left").dropna(subset=["group"])
    snvs["mag_id"] = snvs["bin"].astype(str)

    if args.analysis == "dist":
        # --- Distribution path ---
        if args.freq_col not in snvs.columns:
            raise SystemExit(f"--analysis=dist requires '{args.freq_col}' in SNVs")
        vtab = build_variant_table(snvs, args.freq_col)
        
        # Apply feature-level filtering (same as GLM path)
        log.info("Applying feature filtering to distribution analysis...")
        by_feat = vtab.groupby("feature")
        keep_mask = (by_feat["mag_id"].nunique() >= args.min_strata) & (by_feat["iso_id"].nunique() >= args.min_isolates)
        keep_feats = keep_mask[keep_mask].index.tolist()
        log.info(f"Distribution features before filtering: {len(by_feat)}")
        log.info(f"Distribution features meeting criteria: {len(keep_feats)}")
        
        vtab = vtab[vtab["feature"].isin(keep_feats)].copy()
        
        results = []
        for feat, dfv in vtab.groupby("feature", sort=False):
            info_list = []
            if args.dist_test in ("logit_lm","both"):
                info_list.append(fit_feature_logit_lm(dfv))
            if args.dist_test in ("qreg80","both"):
                info_list.append(fit_feature_qreg(dfv, tau=args.tau))
            for info in info_list:
                if info:
                    info.update({"feature": feat})
                    results.append(info)
        if not results:
            raise SystemExit("No features produced a valid distribution fit.")
        res = pd.DataFrame(results)

        # FDR for each p_* column
        for c in [c for c in res.columns if c.startswith("p_")]:
            res[c.replace("p_","FDR_")] = multipletests(res[c].astype(float), method="fdr_bh")[1]

        out_csv = out_dir / f"ibd_vs_nonibd_{args.family}_dist_{args.dist_test}.csv"
        res.to_csv(out_csv, index=False)
        log.info("✓ Wrote %s (%d features)", out_csv, len(res))
        return  # stop after distribution analysis

    if args.analysis == "qc":
        log.info("Running simple distribution tests (patient-isolate level)")
        # For dist_simple, we need freq_col
        cols_needed = ["patient_id","bin","is_sweep","locus_tag",args.freq_col]
        missing = set(cols_needed) - set(snvs.columns)
        if missing:
            raise ValueError(f"QC analysis requires columns: {missing}")
        
        # Check if group column already exists from earlier processing
        if 'group' in snvs.columns:
            log.info(f"SNVs already has group column with {snvs['group'].nunique()} unique groups")
            log.info(f"Group value counts: {snvs['group'].value_counts()}")
            # Just drop rows without group info
            snvs = snvs.dropna(subset=["group"])
            log.info(f"SNVs with group info: {len(snvs)} rows")
        else:
            # Process metadata to get group information (reload original metadata)
            meta_original = read_table_auto(args.meta_file)
            meta = map_dx_to_group(meta_original, args.grouping)
            log.info(f"Meta with groups: {meta.shape}")
            
            # Make sure group is attached
            snvs = snvs.merge(meta, on="patient_id", how="left").dropna(subset=["group"])
            log.info(f"SNVs with group info: {len(snvs)} rows")
        
        # Build iso×feature presence map
        present_iso_feat = build_iso_feature_presence(mag_genes, args.family)
        log.info(f"Built iso×feature presence map: {len(present_iso_feat)} combinations")
        
        # Build patient-isolate-burden table (presence-restricted)
        burden_df = build_patient_isolate_burden_table(snvs, args.freq_col, present_iso_feat)
        log.info(f"Patient-isolate-burden table built: {len(burden_df)} rows (presence-restricted)")
        log.info(f"Burden table columns: {list(burden_df.columns)}")
        
        # Run burden tests
        res = run_burden_tests(burden_df, 
                              min_patients_per_group=args.min_patients_per_group,
                              min_signal_pct=args.min_signal_pct,
                              q=args.q)
        
        if res.empty:
            raise SystemExit("No features passed filters for dist_simple.")

        log.info(f"Burden test results: {len(res)} rows")
        log.info(f"Metrics tested: {res['metric'].unique()}")
        
        # FDR per comparison and metric
        for comp in res["comparison"].unique():
            for metric in res["metric"].unique():
                mask = (res["comparison"] == comp) & (res["metric"] == metric)
                if mask.sum() > 0:
                    res.loc[mask, "FDR"] = multipletests(res.loc[mask,"p"], method="fdr_bh")[1]

        out_csv = out_dir / f"ibd_vs_nonibd_{args.family}_qc.csv"
        res.to_csv(out_csv, index=False)
        log.info(f"✓ Wrote {out_csv} ({len(res)} rows)")
        return

    # Build per-isolate×feature table
    log.info("Building per-isolate×feature aggregation table...")
    opp_gene_df = read_table_auto(args.opp_genes_file) if args.opp_genes_file else None
    if opp_gene_df is not None:
        log.info(f"Loaded opportunity gene file: {len(opp_gene_df)} rows")
    
    # For GLM path, snvs already has group attached
    agg = build_agg(snvs, meta, args.family, args.opp, opp_gene_df, args.y_metric, args.freq_col)
    if agg.empty: raise SystemExit("No per-feature isolates after preprocessing.")
    log.info(f"Aggregation table built: {len(agg)} rows")
    log.info(f"Aggregation columns: {list(agg.columns)}")

    # Quick diagnostic to understand FE saturation
    diag = agg.groupby("feature").agg(
        n_obs=("iso_id","size"),
        n_mag=("mag_id","nunique"),
        min_iso_per_mag=("iso_id", lambda x: pd.Series(x).groupby(agg.loc[x.index,"mag_id"]).nunique().min())
    ).reset_index()
    log.info("Feature diagnostics (top 10 by obs count):")
    log.info(f"{diag.sort_values('n_obs', ascending=False).head(10).to_string()}")

    # Basic filtering (ensure proper comparison set)
    log.info("Applying feature filtering...")
    by_feat = agg.groupby("feature")
    keep_mask = (by_feat["mag_id"].nunique() >= args.min_strata) & (by_feat["iso_id"].nunique() >= args.min_isolates)
    keep_feats = keep_mask[keep_mask].index.tolist()
    log.info(f"Features before filtering: {len(by_feat)}")
    log.info(f"Features meeting criteria: {len(keep_feats)}")
    
    agg = agg[agg["feature"].isin(keep_feats)].copy()
    log.info("Features kept after min-strata=%d & min-isolates=%d: %d",
             args.min_strata, args.min_isolates, len(keep_feats))

    # Fit per-feature
    log.info(f"Starting GLM fitting for {len(keep_feats)} features...")
    results = []
    for i, (feat, dko) in enumerate(agg.groupby("feature", sort=False)):
        if i % 10 == 0:
            log.info(f"Fitting feature {i+1}/{len(keep_feats)}: {feat}")
        info = fit_one_feature(dko, engine=args.engine, sparse_cutoff=args.sparse_cutoff, 
                              y_metric=args.y_metric, robust_se=args.robust_se, use_gee=args.gee_engine,
                              tweedie_p=args.tweedie_p)
        if info is None: 
            log.debug(f"Feature {feat} failed to fit")
            continue
        info.update({"feature": feat})
        results.append(info)

    if not results:
        raise SystemExit("No features produced a valid fit. Try lowering filters.")

    log.info(f"Successfully fitted {len(results)} features")
    res = pd.DataFrame(results)
    log.info(f"Results DataFrame shape: {res.shape}")
    log.info(f"Results columns: {list(res.columns)}")
    
    # multiple testing
    log.info("Applying multiple testing correction...")
    
    # Handle both groupwise (freqsum + multi-group hits) and single-group outputs
    pcols = [c for c in res.columns if c.startswith("p_")]
    if pcols:
        # Groupwise output (freqsum or multi-group hits)
        for c in pcols:
            res[c.replace("p_","FDR_")] = multipletests(res[c].astype(float), method="fdr_bh")[1]
        # Create log2RR columns for each group
        for c in [c for c in res.columns if c.startswith("RR_")]:
            lev = c.replace("RR_", "")
            res[f"log2RR_{lev}"] = np.log2(res[c])
    else:
        # Single-group output (legacy hits model)
        res["p"] = res["p"].astype(float)
        res["FDR"] = multipletests(res["p"], method="fdr_bh")[1]
        res["log2RR"] = np.log2(res["RR"])
    
    # Check significance
    if pcols:
        # Groupwise output
        sig_counts = {}
        for c in pcols:
            lev = c.replace("p_", "")
            sig_counts[lev] = {
                "p<0.05": (res[c] < 0.05).sum(),
                "p<0.01": (res[c] < 0.01).sum(),
                "FDR<0.05": (res[c.replace("p_","FDR_")] < 0.05).sum()
            }
        log.info(f"Significant features by group: {sig_counts}")
    else:
        # Single-group output
        sig_05 = (res["p"] < 0.05).sum()
        sig_01 = (res["p"] < 0.01).sum()
        fdr_05 = (res["FDR"] < 0.05).sum()
        log.info(f"Significant features: p<0.05: {sig_05}, p<0.01: {sig_01}, FDR<0.05: {fdr_05}")
    
    # tidy cols - handle both single and multi-group cases
    groupwise = any(c.startswith("RR_") for c in res.columns)
    
    if groupwise:
        base_cols = ["feature", "effect", "n_iso", "n_mag"]
        if "total_sweeps" in res.columns: base_cols.append("total_sweeps")
        if "total_freqsum" in res.columns: base_cols.append("total_freqsum")
        for pref in ["RR","p","FDR","ci_low","ci_high","log2RR"]:
            base_cols += [c for c in res.columns if c.startswith(pref + "_")]
    else:
        base_cols = ["feature","effect","RR","log2RR","p","FDR","ci_low","ci_high",
                     "n_iso","n_mag","total_sweeps"]
        if "total_freqsum" in res.columns:
            base_cols.append("total_freqsum")
            
    # Only keep columns that exist
    base_cols = [c for c in base_cols if c in res.columns]
    res = res[base_cols]
    
    # Sort by FDR if it exists, otherwise by first p column
    if "FDR" in res.columns:
        res = res.sort_values("FDR")
    elif any(c.startswith("FDR_") for c in res.columns):
        fdr_cols = [c for c in res.columns if c.startswith("FDR_")]
        res = res.sort_values(fdr_cols[0])

    # write
    out_csv = out_dir / f"ibd_vs_nonibd_{args.family}.csv"
    res.to_csv(out_csv, index=False)
    log.info("✓ Wrote %s (%d features)", out_csv, len(res))

    with open(out_dir / "summary.json", "w") as fh:
        json.dump({
            "family": args.family,
            "engine": args.engine,
            "opp": args.opp,
            "parquet_dir": str(args.parquet_dir),
            "n_features_fit": int(len(res)),
            "n_isolates": int(agg["iso_id"].nunique()),
            "n_mags": int(agg["mag_id"].nunique()),
            "filters": {
                "min_strata": args.min_strata,
                "min_isolates": args.min_isolates,
                "sparse_cutoff": args.sparse_cutoff
            }
        }, fh, indent=2)

    if args.volcano:
        # For volcano plots, use the first available log2RR and p columns
        if "log2RR" in res.columns and "p" in res.columns:
            volcano(res, out_dir / f"volcano_{args.family}.png",
                    f"{args.family.upper()} – IBD vs nonIBD")
        else:
            log.info("Skipping volcano plot - no single log2RR/p columns available for multi-group analysis")

if __name__ == "__main__":
    main()
