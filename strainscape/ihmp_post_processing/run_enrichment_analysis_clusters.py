#!/usr/bin/env python3
"""
Feature enrichment analysis for swept genes with length bias correction.

Biological question: Which functional categories are enriched among genes 
experiencing selective sweeps (high-frequency variants)?

Statistical approach:
  - Binary outcome: gene swept (≥1 SNV with freq_range ≥ freq_cut) vs not swept
  - Null: Genes are randomly selected with probability proportional to gene length
  - Test: Modified hypergeometric test with length-bias correction (GOseq-style)
  
Key innovations vs original:
  1. Uses only sweep events (not all SNVs)
  2. Binary gene-level analysis (avoids within-gene linkage)
  3. Gene length correction via probability weighting function (PWF)
  4. Proper background = all genes in MAG with features

References:
  - Young et al. (2010) Genome Biology - GOseq method
  - Mi et al. (2012) PLOS ONE - GOglm logistic regression approach

Input files:
  --snv-path      Parquet/CSV file or directory; must include: patient_id, bin, freq_range, locus_tag
  --meta-file     CSV/TSV with patient_id, diagnosis (CD/UC/nonIBD) to derive group=IBD/nonIBD
  --mag-genes     TSV compiled earlier with columns:
                   patient_id  bin  scaffold  locus_tag  cluster_rep  rep_gene  rep_dbxrefs
  --filter-file   (optional) TSV with patient_id  bin     → restrict analysis to these isolates

Output:
  sweep_enrichment_{feature_type}_{group}.csv  (all features)
  sweep_enrichment_{feature_type}_{group}_significant.csv  (FDR<..., min participants)
  run_params.json

Usage example:
  python strainscape/ihmp_post_processing/run_enrichment_analysis_clusters.py \
    --snv-path /Users/reneeoles/Desktop/strainscape_output/iHMP/output/all_snvs \
    --meta-file  /Users/reneeoles/Desktop/strainscape_output/iHMP/metadata/meta_map.csv \
    --mag-genes /Users/reneeoles/Desktop/strainscape_output/iHMP/input/all_genes.tsv \
    --out-dir /Users/reneeoles/Desktop/strainscape_output/iHMP/output/enrichment_analysis \
    --filter-file /Users/reneeoles/Desktop/strainscape_output/iHMP/output/figures/sweeps/patient_bin_1000.txt
"""

from __future__ import annotations
import argparse, json, logging, re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import UnivariateSpline
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

log = logging.getLogger("sweep_enrichment")
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# ─────────────────────────────────────────────────────────────────────────────
# Parsers for features coming from MAG gene file
# ─────────────────────────────────────────────────────────────────────────────

KO_PAT = re.compile(r'^K\d{5}$')
GO_PAT = re.compile(r'^\d{7}$')
EC_PAT = re.compile(r'^\d+\.\d+\.\d+\.\d+$')

MGE_GENES = {
        # Conjugative transfer genes
        'traA', 'traB', 'traC', 'traD', 'traE', 'traF', 'traG', 'traH', 'traI', 'traJ',
        'traK', 'traL', 'traM', 'traN', 'traO', 'traP', 'traQ', 'traR', 'traS', 'traT',
        'traU', 'traV', 'traW', 'traX', 'traY', 'traZ',
        # Mobile element mobilization
        'mobA', 'mobB', 'mobC', 'mobD', 'mobE', 'mobF', 'mobG', 'mobH', 'mobI', 'mobJ',
        'mobK', 'mobL', 'mobM', 'mobN', 'mobO', 'mobP', 'mobQ', 'mobR', 'mobS', 'mobT',
        'mobU', 'mobV', 'mobW', 'mobX', 'mobY', 'mobZ',
        # Site-specific recombination
        'xerA', 'xerB', 'xerC', 'xerD', 'xerE', 'xerF', 'xerG', 'xerH', 'xerI', 'xerJ',
        'xerK', 'xerL', 'xerM', 'xerN', 'xerO', 'xerP', 'xerQ', 'xerR', 'xerS', 'xerT',
        'virD4', 'virB',
        # Transposons and IS elements
        'tnpA', 'tnpB', 'tnpC', 'tnpD', 'tnpE', 'tnpF', 'tnpG', 'tnpH', 'tnpI', 'tnpJ',
        'tnpK', 'tnpL', 'tnpM', 'tnpN', 'tnpO', 'tnpP', 'tnpQ', 'tnpR', 'tnpS', 'tnpT',
        'tnpU', 'tnpV', 'tnpW', 'tnpX', 'tnpY', 'tnpZ',
        'is1', 'is2', 'is3', 'is4', 'is5', 'is6', 'is7', 'is8', 'is9', 'is10',
        'is11', 'is12', 'is13', 'is14', 'is15', 'is16', 'is17', 'is18', 'is19', 'is20',
        'is21', 'is22', 'is23', 'is24', 'is25', 'is26', 'is27', 'is28', 'is29', 'is30',
        # Integrases and recombinases
        'int', 'intA', 'intB', 'intC', 'intD', 'intE', 'intF', 'intG', 'intH', 'intI',
        'intJ', 'intK', 'intL', 'intM', 'intN', 'intO', 'intP', 'intQ', 'intR', 'intS',
        'intT', 'intU', 'intV', 'intW', 'intX', 'intY', 'intZ',
        'recA', 'recB', 'recC', 'recD', 'recE', 'recF', 'recG', 'recH', 'recI', 'recJ',
        'recK', 'recL', 'recM', 'recN', 'recO', 'recP', 'recQ', 'recR', 'recS', 'recT',
        # Phage-related genes
        'phage', 'phage_integrase', 'phage_repressor', 'phage_tail', 'phage_head',
        'phage_portal', 'phage_terminase', 'phage_major_capsid', 'phage_minor_capsid',
        # Plasmid maintenance
        'parA', 'parB', 'parC', 'parD', 'parE', 'parF', 'parG', 'parH', 'parI', 'parJ',
        'parK', 'parL', 'parM', 'parN', 'parO', 'parP', 'parQ', 'parR', 'parS', 'parT',
        # Replication initiation
        'repA', 'repB', 'repC', 'repD', 'repE', 'repF', 'repG', 'repH', 'repI', 'repJ',
        'repK', 'repL', 'repM', 'repN', 'repO', 'repP', 'repQ', 'repR', 'repS', 'repT',
        # Partitioning and stability
        'stbA', 'stbB', 'stbC', 'stbD', 'stbE', 'stbF', 'stbG', 'stbH', 'stbI', 'stbJ',
        'stbK', 'stbL', 'stbM', 'stbN', 'stbO', 'stbP', 'stbQ', 'stbR', 'stbS', 'stbT',
        # Resolution systems
        'resA', 'resB', 'resC', 'resD', 'resE', 'resF', 'resG', 'resH', 'resI', 'resJ',
        'resK', 'resL', 'resM', 'resN', 'resO', 'resP', 'resQ', 'resR', 'resS', 'resT',
        'nan', 'rteC',
        # Ribosomal genes
        'virD4', 'rpoA', 'rpoB', 'rpoC', 'rpoD', 'rpoE', 'rpoF', 'rpoG', 'rpoH', 'rpoI', 'rpoJ',
        'rpoK', 'rpoL', 'rpoM', 'rpoN', 'rpoO', 'rpoP', 'rpoQ', 'rpoR', 'rpoS', 'rpoT',
        'rpoU', 'rpoV', 'rpoW', 'rpoX', 'rpoY', 'rpoZ', 'rplV', 'rpsS', 'rpsB', 'rpsA',
        'prsM', 'rpsS', 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ',
        'rplC', 'rpsM', 'rplL', 'rplK', 'rplJ', 'rplI', 'rplH', 'rplG', 'rplF', 'rplE', 
        'rplD', 'rplC', 'rplB', 'rplA', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ', 'rplR',
        'rplS', 'rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY', 'rplZ'
    }

def _split_dbx(x):
    if pd.isna(x) or x is None: return []
    return [t for t in re.split(r"[,;|\s]+", str(x).strip()) if t and t.lower() not in {"nan","none"}]


def feats_from_dbxrefs(dbx, family):
    vals=[]
    for t in _split_dbx(dbx):
        u=t.upper()
        if family=="kegg":
            u = u.split(":",1)[1] if u.startswith("KEGG:") else u
            if KO_PAT.match(u): vals.append("kegg:"+u)
        elif family=="go":
            u = u.split(":",1)[1] if u.startswith("GO:") else u
            if GO_PAT.match(u): vals.append("go:"+u)
        elif family=="ec":
            u = u.split(":",1)[1] if u.startswith("EC:") else u
            if EC_PAT.match(u): vals.append("ec:"+u)
    return list(dict.fromkeys(vals))

def feats_from_gene(rep_gene):
    if pd.isna(rep_gene) or not str(rep_gene).strip(): return []
    toks = [t.strip() for t in re.split(r"[ ,;|]+", str(rep_gene)) if t.strip()]
    bad = {gene.lower() for gene in MGE_GENES}
    return [t for t in toks if t.lower() not in bad]

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_meta(meta_path: Path) -> pd.DataFrame:
    log.info(f"Loading metadata from: {meta_path}")
    try:
        meta = pd.read_csv(meta_path, sep=None, engine="python")
        log.info(f"Loaded metadata with auto-detected separator. Shape: {meta.shape}")
    except Exception:
        meta = pd.read_csv(meta_path, sep="\t")
        log.info(f"Loaded metadata with tab separator. Shape: {meta.shape}")
    
    log.info(f"Metadata columns: {list(meta.columns)}")
    
    # harmonize ids + groups
    cols = {c.lower(): c for c in meta.columns}
    if "patient_id" not in cols:
        if "participant id" in cols: 
            meta = meta.rename(columns={cols["participant id"]:"patient_id"})
            log.info("Renamed 'Participant ID' to 'patient_id'")
        elif "run" in cols:         
            meta = meta.rename(columns={cols["run"]:"patient_id"})
            log.info("Renamed 'Run' to 'patient_id'")
    if "diagnosis" not in meta.columns:
        raise ValueError("metadata needs a 'diagnosis' column")
    
    meta["patient_id"] = meta["patient_id"].astype(str).str.strip()
    dx = meta["diagnosis"].astype(str).str.strip().str.lower()
    log.info(f"Diagnosis distribution: {dx.value_counts().to_dict()}")
    
    mapping = {
        "cd":"IBD","crohn":"IBD","crohn's":"IBD","crohns":"IBD",
        "uc":"IBD","ulcerative colitis":"IBD",
        "nonibd":"nonIBD","non-ibd":"nonIBD","healthy":"nonIBD","control":"nonIBD"
    }
    meta["group"] = dx.map(mapping)
    log.info(f"Group distribution after mapping: {meta['group'].value_counts().to_dict()}")
    
    meta = meta[["patient_id","group"]].dropna().drop_duplicates()
    log.info(f"Final metadata shape after cleanup: {meta.shape}")
    return meta

def read_table_auto(p: Path) -> pd.DataFrame:
    with open(p, 'r') as f: sep = ("\t" if "\t" in f.readline() else ",")
    return pd.read_csv(p, sep=sep)

def load_snvs(snv_path: Path) -> pd.DataFrame:
    cols = ["patient_id","bin","locus_tag","freq_range","type","start","stop"]
    if snv_path.is_file():
        df = pd.read_parquet(snv_path) if snv_path.suffix.lower()==".parquet" else read_table_auto(snv_path)
    else:
        files = sorted(list(snv_path.glob("*.parquet"))) or sorted(list(snv_path.glob("*.csv")))
        df = pd.concat([pd.read_parquet(p) if p.suffix==".parquet" else read_table_auto(p) for p in files], ignore_index=True)
    
    # Filter to only CDS (protein-coding) genes
    if "type" in df.columns:
        log.info("Filtering to CDS (protein-coding) genes only...")
        before_count = len(df)
        df = df[df["type"] == "cds"].copy()
        after_count = len(df)
        log.info(f"CDS filtering: {before_count} -> {after_count} SNVs ({after_count/before_count:.1%} retained)")
    else:
        log.warning("No 'type' column found - skipping CDS filtering")
    
    # Filter out mobile genetic elements (MGEs) if a gene name column is present
    if 'gene' in df.columns:
        log.info("Filtering out mobile genetic elements...")
        mge_genes_lower = {gene.lower() for gene in MGE_GENES}
        before_mge_count = len(df)
        df = df[~df['gene'].str.lower().isin(mge_genes_lower)].copy()
        after_mge_count = len(df)
        log.info(f"MGE filtering: {before_mge_count} -> {after_mge_count} SNVs ({after_mge_count/before_mge_count:.1%} retained)")
    else:
        log.info("No 'gene' column found - skipping MGE filtering")
    
    
    # normalize
    df["patient_id"] = df["patient_id"].astype(str).str.strip()
    df["bin"] = df["bin"].astype(str).str.strip()
    df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
    df["freq_range"] = pd.to_numeric(df["freq_range"], errors="coerce").fillna(0.0).clip(0,1)
    df["iso_id"] = df["patient_id"] + "|" + df["bin"]
    return df[["patient_id","bin","iso_id","locus_tag","freq_range","start","stop"]]

def split_locus_tags(df: pd.DataFrame, col="locus_tag", weight_col="freq_range", sep="/") -> pd.DataFrame:
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

def load_mag_genes(mag_path: Path) -> pd.DataFrame:
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

# ─────────────────────────────────────────────────────────────────────────────
# Core functions - Length-bias corrected enrichment analysis
# ─────────────────────────────────────────────────────────────────────────────
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

def calculate_gene_lengths(mag_genes: pd.DataFrame, snvs: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate gene length from SNV data coordinates.
    Falls back to default length if coordinates unavailable.
    """
    mg = mag_genes.copy()
    
    # Extract gene lengths from SNV data
    log.info("Extracting gene lengths from SNV data...")
    gene_lengths = extract_gene_lengths_from_snvs(snvs)
    log.info(f"Extracted lengths for {len(gene_lengths)} genes from SNV data")
    
    # Merge gene lengths with mag_genes
    mg = mg.merge(gene_lengths, on="locus_tag", how="left")
    
    # Fill missing lengths with default
    default_gene_len = 1000.0
    mg["gene_len"] = mg["gene_len"].fillna(default_gene_len)
    
    log.info(f"Gene length statistics: mean={mg['gene_len'].mean():.1f}, median={mg['gene_len'].median():.1f}, range={mg['gene_len'].min():.0f}-{mg['gene_len'].max():.0f}")
    
    # Handle invalid lengths
    mg["gene_len"] = mg["gene_len"].clip(lower=100)
    
    return mg

def fit_length_bias_pwf(swept_genes: pd.DataFrame, 
                        all_genes: pd.DataFrame,
                        method: str = 'logistic') -> pd.DataFrame:
    """
    Fit probability weighting function (PWF) to model relationship between
    gene length and probability of being swept.
    
    This is the key to accounting for length bias.
    
    Args:
        swept_genes: genes with ≥1 sweep (columns: iso_id, locus_tag)
        all_genes: all genes with features (columns: iso_id, locus_tag, gene_len)
        method: 'logistic' (recommended) or 'spline' (GOseq-style)
    
    Returns:
        DataFrame with columns: iso_id, locus_tag, gene_len, is_swept, pwf
    """
    # Merge to get swept status
    genes = all_genes.copy()
    genes['is_swept'] = genes.apply(
        lambda r: int((r['iso_id'], r['locus_tag']) in 
                     set(zip(swept_genes['iso_id'], swept_genes['locus_tag']))),
        axis=1
    )
    
    log.info(f"PWF fitting: {genes['is_swept'].sum()} swept / {len(genes)} total genes")
    
    if method == 'logistic':
        # Logistic regression: logit(P(swept)) ~ log(length)
        # Following Mi et al. 2012 GOglm approach
        from sklearn.linear_model import LogisticRegression
        
        X = np.log(genes['gene_len'].values).reshape(-1, 1)
        y = genes['is_swept'].values
        
        # Fit with balanced class weights
        model = LogisticRegression(class_weight='balanced', max_iter=1000)
        model.fit(X, y)
        
        # Predict probabilities
        pwf = model.predict_proba(X)[:, 1]
        genes['pwf'] = pwf
        
        log.info(f"Logistic PWF: coef={model.coef_[0][0]:.4f}, "
                f"mean_pwf={pwf.mean():.4f}")
        
    elif method == 'spline':
        # Monotonic spline following Young et al. 2010 GOseq
        # Bin genes by length and calculate proportion swept per bin
        genes = genes.sort_values('gene_len')
        genes['length_bin'] = pd.qcut(genes['gene_len'], q=20, duplicates='drop')
        
        bin_stats = genes.groupby('length_bin').agg({
            'is_swept': 'mean',
            'gene_len': 'median'
        }).reset_index()
        
        # Fit monotonic spline
        spline = UnivariateSpline(
            np.log(bin_stats['gene_len']),
            bin_stats['is_swept'],
            s=0.5,  # smoothing
            k=3     # cubic
        )
        
        genes['pwf'] = spline(np.log(genes['gene_len']))
        genes['pwf'] = genes['pwf'].clip(0.001, 0.999)  # bound away from 0/1
        
    else:
        raise ValueError(f"Unknown PWF method: {method}")
    
    return genes[['iso_id', 'locus_tag', 'gene_len', 'is_swept', 'pwf']]

def wallenius_test(k: int, n: int, m: int, M: int, 
                   weights: np.ndarray, odds: float = 1.0) -> float:
    """
    Wallenius non-central hypergeometric test for enrichment.
    
    Args:
        k: # swept genes in feature
        n: # total genes in feature
        m: # swept genes genome-wide
        M: # total genes genome-wide
        weights: length-based weights for all genes
        odds: odds ratio (typically 1 for null)
    
    Returns:
        p-value for one-sided enrichment test
    """
    # For large samples, use normal approximation
    # E[X] and Var[X] under Wallenius non-central hypergeometric
    
    # Expected value
    p_select = m / M
    E_X = n * p_select
    
    # Variance (approximation)
    Var_X = n * p_select * (1 - p_select) * (M - n) / (M - 1)
    
    # Z-score
    z = (k - E_X) / np.sqrt(Var_X + 1e-10)
    
    # One-sided p-value
    p = 1 - stats.norm.cdf(z)
    
    return max(p, 1e-300)

def feature_enrichment_lengthcorrected(
    snvs: pd.DataFrame,
    mag_genes: pd.DataFrame,
    meta: pd.DataFrame,
    family: str,
    freq_cut: float = 0.6,
    min_swept_genes: int = 10,
    min_participants: int = 10
) -> pd.DataFrame:
    """
    Test for feature enrichment among swept genes with length bias correction.
    
    Returns:
        DataFrame with columns: feature, k, n, m, M, OR, p, FDR, n_participants
        where:
          k = # swept genes with this feature
          n = # total genes with this feature
          m = # swept genes genome-wide
          M = # total genes genome-wide
    """
    # 1) Identify swept genes (binary: gene has ≥1 high-freq SNV)
    snvs['iso_id'] = snvs['patient_id'].astype(str) + '|' + snvs['bin'].astype(str)
    snvs['is_high_freq'] = snvs['freq_range'] >= freq_cut
    
    swept_genes = (snvs[snvs['is_high_freq']]
                   .groupby(['iso_id', 'locus_tag'])
                   .size()
                   .reset_index(name='n_sweeps'))
    swept_genes = swept_genes[swept_genes['n_sweeps'] > 0]
    
    log.info(f"Identified {len(swept_genes)} swept genes across {swept_genes['iso_id'].nunique()} isolates")
    
    # 2) Calculate gene lengths from SNV data
    mag_genes = calculate_gene_lengths(mag_genes, snvs)
    
    # 3) Build feature map (gene -> features)
    if family == 'gene':
        feat_series = mag_genes['rep_gene'].apply(
            lambda x: feats_from_gene(x) if pd.notna(x) else []
        )
    else:
        feat_series = mag_genes['rep_dbxrefs'].apply(
            lambda x: feats_from_dbxrefs(x, family) if pd.notna(x) else []
        )
    
    fm = (mag_genes.assign(feature=feat_series)
          .explode('feature')
          .dropna(subset=['feature'])
          .loc[:, ['iso_id', 'locus_tag', 'feature', 'gene_len']])
    
    # 4) Fit PWF to account for length bias
    all_genes_with_features = fm[['iso_id', 'locus_tag', 'gene_len']].drop_duplicates()
    pwf_data = fit_length_bias_pwf(swept_genes, all_genes_with_features, method='logistic')
    
    # 5) Merge swept status with features
    fm = fm.merge(pwf_data[['iso_id', 'locus_tag', 'is_swept', 'pwf']], 
                  on=['iso_id', 'locus_tag'], how='left')
    
    # 6) Per-feature statistics
    results = []
    for feat, grp in fm.groupby('feature'):
        # Counts
        k = grp['is_swept'].sum()  # swept genes with this feature
        n = len(grp)               # total genes with this feature
        
        if k < min_swept_genes:
            continue
        
        # Genome-wide totals (from pwf_data, not grp, to avoid biased denominator)
        m = pwf_data['is_swept'].sum()  # total swept genes
        M = len(pwf_data)                # total genes
        
        # Participant support
        n_participants = (grp[grp['is_swept'] == 1]
                         .assign(patient_id=lambda x: x['iso_id'].str.split('|').str[0])
                         ['patient_id'].nunique())
        
        if n_participants < min_participants:
            continue
        
        # Length-corrected test using Wallenius approximation
        weights = pwf_data['pwf'].values
        p = wallenius_test(k, n, m, M, weights, odds=1.0)
        
        # Odds ratio: observed/expected
        expected = n * (m / M)
        OR = (k + 0.5) / (expected + 0.5)
        
        results.append({
            'feature': feat,
            'k_swept': int(k),
            'n_total': int(n),
            'm_swept_genome': int(m),
            'M_total_genome': int(M),
            'OR': float(OR),
            'p_value': float(p),
            'n_participants': int(n_participants)
        })
    
    if not results:
        return pd.DataFrame()
    
    res = pd.DataFrame(results)
    
    # FDR correction
    res['FDR'] = multipletests(res['p_value'], method='fdr_bh')[1]
    
    # Sort by significance
    res = res.sort_values(['FDR', 'p_value'])
    
    log.info(f"Tested {len(res)} features, {(res['FDR'] < 0.05).sum()} significant at FDR<0.05")
    
    return res

# Helper functions (same as original)
def feats_from_gene(rep_gene):
    if pd.isna(rep_gene) or not str(rep_gene).strip(): return []
    toks = [t.strip() for t in re.split(r"[ ,;|]+", str(rep_gene)) if t.strip()]
    bad = {"hypothetical","hypothetical_protein","uncharacterized","putative","unknown","none","nan"}
    return [t for t in toks if t.lower() not in bad]

def feats_from_dbxrefs(dbx, family):
    if pd.isna(dbx): return []
    vals = []
    for t in re.split(r"[,;|\s]+", str(dbx).strip()):
        u = t.upper()
        if family == "kegg":
            u = u.split(":",1)[1] if u.startswith("KEGG:") else u
            if KO_PAT.match(u): vals.append(u)
        elif family == "go":
            u = u.split(":",1)[1] if u.startswith("GO:") else u
            if GO_PAT.match(u): vals.append(u)
        elif family == "ec":
            u = u.split(":",1)[1] if u.startswith("EC:") else u
            if EC_PAT.match(u): vals.append(u)
    return list(dict.fromkeys(vals))

# ─────────────────────────────────────────────────────────────────────────────
# Core builders - Legacy (for reference)
# ─────────────────────────────────────────────────────────────────────────────

def build_feature_map(mag_genes: pd.DataFrame, family: str) -> pd.DataFrame:
    mg = mag_genes.copy()
    mg["iso_id"] = mg["patient_id"].astype(str)+"|"+mg["bin"].astype(str)
    if family=="gene":
        feats = mg["rep_gene"].apply(feats_from_gene)
    else:
        feats = mg["rep_dbxrefs"].apply(lambda s: feats_from_dbxrefs(s, family))
    fm = (mg.assign(feature_list=feats)
            .explode("feature_list")
            .dropna(subset=["feature_list"]))
    fm = fm[["iso_id","locus_tag","feature_list"]].rename(columns={"feature_list":"feature"})
    # per-gene degree & split weight
    deg = fm.groupby(["iso_id","locus_tag"], observed=True)["feature"].nunique().rename("deg").reset_index()
    fm = fm.merge(deg, on=["iso_id","locus_tag"], how="left")
    fm["w_split"] = 1.0 / fm["deg"].clip(lower=1)
    return fm  # columns: iso_id, locus_tag, feature, w_split

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run(args):
    log.info("="*60)
    log.info("STARTING SWEEP ENRICHMENT ANALYSIS WITH LENGTH BIAS CORRECTION")
    log.info("="*60)
    log.info(f"Output directory: {args.out_dir}")
    log.info(f"Feature type: {args.feature_type}")
    log.info(f"Group filter: {args.group}")
    log.info(f"Frequency cutoff: {args.freq_cut}")
    
    out_dir = args.out_dir; out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading input data...")
    meta = load_meta(args.meta_file)
    snv  = load_snvs(args.snv_path)
    mg   = load_mag_genes(args.mag_genes)
    flt  = load_filter(args.filter_file)
    
    # CRITICAL: Split locus tags containing '/' into separate rows and divide freq_range
    log.info("Splitting locus tags and dividing freq_range weights...")
    snv = split_locus_tags(snv, weight_col="freq_range")

    # Apply filtering if provided
    if flt is not None:
        log.info("Applying patient+bin filter...")
        log.info(f"Filter contains {len(flt)} iso_ids")
        log.info(f"SNVs before filtering: {len(snv)} rows")
        log.info(f"MAG genes before filtering: {len(mg)} rows")
        snv = snv[snv["iso_id"].isin(flt["iso_id"])]
        mg  = mg[mg["iso_id"].isin(flt["iso_id"])]
        log.info(f"SNVs after filtering: {len(snv)} rows")
        log.info(f"MAG genes after filtering: {len(mg)} rows")
        log.info(f"Unique iso_ids in filtered SNVs: {snv['iso_id'].nunique()}")
        log.info(f"Unique iso_ids in filtered MAG genes: {mg['iso_id'].nunique()}")
    else:
        log.info("No filter provided - using all data")

    # Apply group filtering
    if args.group != "ALL":
        log.info(f"Filtering to group: {args.group}")
        snv = snv.merge(meta, on="patient_id", how="left").dropna(subset=["group"])
        snv = snv[snv["group"] == args.group]
        log.info(f"SNVs after group filtering: {len(snv)} rows")

    # Performance optimization: filter MAG genes early to isolates present in SNVs
    log.info("Filtering MAG genes to isolates present in SNVs...")
    keep_iso = snv["iso_id"].unique()
    mg = mg[mg["iso_id"].isin(keep_iso)].copy()
    log.info(f"MAG genes after filtering to SNV isolates: {len(mg)} rows")

    # Run length-corrected enrichment analysis
    log.info("Running sweep enrichment analysis with length bias correction...")
    res = feature_enrichment_lengthcorrected(
        snvs=snv,
        mag_genes=mg,
        meta=meta,
        family=args.feature_type,
        freq_cut=args.freq_cut,
        min_swept_genes=args.min_swept_genes,
        min_participants=args.min_participants
    )

    if res.empty:
        log.warning("No features passed filtering criteria!")
        return

    # Filter significant results
    sig = res[(res["FDR"] < args.p_thresh) &
              (res["OR"] >= args.min_enrichment) &
              (res["n_participants"] >= args.min_participants)].copy()

    stem = f"sweep_enrichment_{args.feature_type}_{args.group.lower()}"
    out_file = out_dir / f"{stem}.csv"
    sig_file = out_dir / f"{stem}_significant.csv"
    
    res.to_csv(out_file, index=False)
    sig.to_csv(sig_file, index=False)
    
    log.info("Wrote: %s (%d features); significant: %s (%d)",
             out_file, len(res), sig_file, len(sig))

    # Between-group comparison (IBD vs non-IBD)
    if args.group == "ALL":
        log.info("Performing between-group comparison (IBD vs non-IBD)...")
        try:
            # Load full data for between-group analysis
            snv_all = load_snvs(args.snv_path)
            snv_all = split_locus_tags(snv_all, weight_col="freq_range")
            
            # Apply same filtering as main analysis
            if flt is not None:
                snv_all = snv_all[snv_all["iso_id"].isin(flt["iso_id"])]
            
            # Run separate analyses for each group
            ibd_snv = snv_all.merge(meta, on="patient_id", how="left")
            ibd_snv = ibd_snv[ibd_snv["group"] == "IBD"]
            nonibd_snv = snv_all.merge(meta, on="patient_id", how="left")
            nonibd_snv = nonibd_snv[nonibd_snv["group"] == "nonIBD"]
            
            ibd_res = feature_enrichment_lengthcorrected(
                snvs=ibd_snv, mag_genes=mg, meta=meta, family=args.feature_type,
                freq_cut=args.freq_cut, min_swept_genes=args.min_swept_genes,
                min_participants=args.min_participants
            )
            nonibd_res = feature_enrichment_lengthcorrected(
                snvs=nonibd_snv, mag_genes=mg, meta=meta, family=args.feature_type,
                freq_cut=args.freq_cut, min_swept_genes=args.min_swept_genes,
                min_participants=args.min_participants
            )
            
            # Merge results for comparison
            ibd_res = ibd_res.rename(columns={'OR': 'OR_IBD', 'FDR': 'FDR_IBD', 'n_participants': 'n_participants_IBD'})
            nonibd_res = nonibd_res.rename(columns={'OR': 'OR_nonIBD', 'FDR': 'FDR_nonIBD', 'n_participants': 'n_participants_nonIBD'})
            
            bt = pd.merge(ibd_res[['feature', 'OR_IBD', 'FDR_IBD', 'n_participants_IBD']], 
                         nonibd_res[['feature', 'OR_nonIBD', 'FDR_nonIBD', 'n_participants_nonIBD']], 
                         on='feature', how='outer').fillna({'OR_IBD': 1.0, 'OR_nonIBD': 1.0})
            
            # Calculate log2 fold change
            bt['log2FC'] = np.log2(bt['OR_IBD'] / bt['OR_nonIBD'])
            
            # Save between-group results
            bt_file = out_dir / f"between_group_{args.feature_type}_IBD_vs_nonIBD.csv"
            bt.to_csv(bt_file, index=False)
            
            log.info("Between-group comparison complete:")
            log.info("  Total features compared: %d", len(bt))
            log.info("  Wrote: %s", bt_file)
            
            if len(bt) > 0:
                log.info("  Top 5 IBD-enriched features: %s", 
                        bt.nlargest(5, 'log2FC')[['feature', 'log2FC', 'OR_IBD', 'OR_nonIBD']].to_string(index=False))
                log.info("  Top 5 non-IBD-enriched features: %s", 
                        bt.nsmallest(5, 'log2FC')[['feature', 'log2FC', 'OR_IBD', 'OR_nonIBD']].to_string(index=False))
                        
        except Exception as e:
            log.warning("Between-group comparison failed: %s", str(e))
    else:
        log.info("Skipping between-group comparison (not running with --group ALL)")

    with open(out_dir / "run_params.json","w") as fh:
        json.dump({
            "feature_type": args.feature_type,
            "group": args.group,
            "method": "sweep_enrichment_length_corrected",
            "freq_cut": args.freq_cut,
            "p_thresh": args.p_thresh,
            "min_enrichment": args.min_enrichment,
            "min_participants": args.min_participants,
            "min_swept_genes": args.min_swept_genes
        }, fh, indent=2)

def parse_args():
    ap = argparse.ArgumentParser("Sweep enrichment analysis with length bias correction")
    ap.add_argument("--snv-path", type=Path, required=True,
                    help="Parquet/CSV file or directory with SNVs (patient_id, bin, freq_range, locus_tag)")
    ap.add_argument("--meta-file", type=Path, required=True,
                    help="Metadata with patient_id, diagnosis")
    ap.add_argument("--mag-genes", type=Path, required=True,
                    help="Compiled MAG gene TSV (patient_id, bin, locus_tag, cluster_rep, rep_gene, rep_dbxrefs)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--group", choices=["ALL","IBD","nonIBD"], default="IBD")
    ap.add_argument("--feature-type", choices=["kegg","go","ec","gene"], default="gene")
    ap.add_argument("--filter-file", type=Path, help="Optional TSV with patient_id and bin to keep")
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--min-enrichment", type=float, default=2.0)
    ap.add_argument("--min-participants", type=int, default=10)
    ap.add_argument("--min-swept-genes", type=int, default=10,
                    help="Minimum # of swept genes per feature to test")
    ap.add_argument("--freq-cut", type=float, default=0.6,
                    help="Frequency cutoff for defining swept genes")
    ap.add_argument("--freq-floor", type=float, default=0.2,
                    help="Minimum freq_range to include (pre-filter low-frequency sites)")
    return ap.parse_args()

if __name__ == "__main__":
    run(parse_args())
