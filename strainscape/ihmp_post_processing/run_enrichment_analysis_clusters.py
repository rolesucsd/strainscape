#!/usr/bin/env python3
"""
Weighted SNV feature enrichment using analytical null expectation (fast, paper-style).

Observed:
  • Use ALL SNVs (not just sweeps), weighted by freq_range.
  • Split multi-locus SNVs and divide freq_range equally across loci.
  • Map each SNV locus -> features via MAG file (rep_gene, rep_dbxrefs).
  • Each SNV contributes its freq_range mass split equally across locus features.
  • Sum across all SNVs to get observed feature mass vector.

Null (analytical):
  • Randomly assign each SNV's mass to a random gene in that MAG (uniform over genes with features).
  • Split by that gene's features using per-gene split weights.
  • Compute expected mass μ and variance analytically using p_i,k and q_i,k.
  • One-sided test for excess mass with FDR correction.

Input files:
  --snv-path      Parquet/CSV file or directory; must include: patient_id, bin, freq_range, locus_tag
  --meta-file     CSV/TSV with patient_id, diagnosis (CD/UC/nonIBD) to derive group=IBD/nonIBD
  --mag-genes     TSV compiled earlier with columns:
                   patient_id  bin  scaffold  locus_tag  cluster_rep  rep_gene  rep_dbxrefs
  --filter-file   (optional) TSV with patient_id  bin     → restrict analysis to these isolates

Output:
  all_snv_weighted_enrichment_{feature_type}_{group}.csv  (all features)
  all_snv_weighted_enrichment_{feature_type}_{group}_significant.csv  (FDR<..., enr≥..., min participants)
  run_params.json

Usage example:
  python run_enrichment_analysis_clusters.py \
    --snv-path /path/snv_parquets \
    --meta-file meta.tsv \
    --mag-genes /path/all_genes.tsv \
    --out-dir out_kegg \
    --feature-type kegg \
    --group IBD \
    --min-participants 3 --min-enrichment 2 --p-thresh 0.05
"""

from __future__ import annotations
import argparse, json, logging, re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

log = logging.getLogger("convergent_mag")
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# ─────────────────────────────────────────────────────────────────────────────
# Parsers for features coming from MAG gene file
# ─────────────────────────────────────────────────────────────────────────────

KO_PAT = re.compile(r'^K\d{5}$')
GO_PAT = re.compile(r'^\d{7}$')
EC_PAT = re.compile(r'^\d+\.\d+\.\d+\.\d+$')

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
    bad = {"hypothetical","hypothetical_protein","uncharacterized","putative","unknown","none","nan"}
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
    cols = ["patient_id","bin","locus_tag","freq_range"]
    if snv_path.is_file():
        df = pd.read_parquet(snv_path) if snv_path.suffix.lower()==".parquet" else read_table_auto(snv_path)
    else:
        files = sorted(list(snv_path.glob("*.parquet"))) or sorted(list(snv_path.glob("*.csv")))
        df = pd.concat([pd.read_parquet(p) if p.suffix==".parquet" else read_table_auto(p) for p in files], ignore_index=True)
    # normalize
    df["patient_id"] = df["patient_id"].astype(str).str.strip()
    df["bin"] = df["bin"].astype(str).str.strip()
    df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
    df["freq_range"] = pd.to_numeric(df["freq_range"], errors="coerce").fillna(0.0).clip(0,1)
    df["iso_id"] = df["patient_id"] + "|" + df["bin"]
    return df[["patient_id","bin","iso_id","locus_tag","freq_range"]]

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
# Core builders - New weighted SNV approach
# ─────────────────────────────────────────────────────────────────────────────

def combine_gene_mass(freqs: pd.Series, how: str = "prob_any", power: float = 1.0) -> float:
    # Emphasize higher-frequency mutations by exponentiating per-site f
    f = np.clip(freqs.astype(float).values, 0.0, 1.0) ** float(power)
    if f.size == 0:
        return 0.0
    if how == "sum":       # sums all site masses (still sensitive to linkage)
        return float(f.sum())
    if how == "max":       # use strongest site as gene mass
        return float(f.max())
    if how == "prob_any":  # 1 - Π(1 - f_j^α): mass of "any mutation in this gene"
        return float(1.0 - np.prod(1.0 - f))
    raise ValueError("how must be one of: prob_any, max, sum")

def observed_feature_mass(snvs: pd.DataFrame, fm: pd.DataFrame,
                          collapse: str = "prob_any",
                          freq_power: float = 1.0
                         ) -> tuple[pd.Series, pd.Series, pd.DataFrame, pd.DataFrame, pd.Series]:
    # Collapse to one mass per (iso_id, locus_tag)
    snvs = snvs.copy()
    gene_mass = (snvs.groupby(["iso_id","locus_tag"], observed=True)["freq_range"]
                      .apply(lambda s: combine_gene_mass(s, how=collapse, power=freq_power))
                      .rename("m_gene").reset_index())

    # Keep only genes that have at least one feature
    joinable = fm[["iso_id","locus_tag"]].drop_duplicates()
    gene_mass = gene_mass.merge(joinable, on=["iso_id","locus_tag"], how="inner", sort=False)

    # Binary "event": gene mutated in this isolate if m_gene > 0
    mutated_genes = gene_mass.loc[gene_mass["m_gene"] > 0, ["iso_id","locus_tag"]].drop_duplicates()

    # Event counts per feature (unique iso×gene events that map to the feature)
    events = (mutated_genes.merge(fm[["iso_id","locus_tag","feature"]], on=["iso_id","locus_tag"], how="inner")
                          .drop_duplicates(["iso_id","locus_tag","feature"])
                          .groupby("feature", observed=True).size().astype(int))

    # Split gene mass across gene's features by w_split
    snv_feat = gene_mass.merge(
        fm[["iso_id","locus_tag","feature","w_split"]],
        on=["iso_id","locus_tag"], how="inner", sort=False
    )
    snv_feat["mass"] = snv_feat["m_gene"] * snv_feat["w_split"]

    # Observed mass per feature
    obs = snv_feat.groupby("feature", observed=True)["mass"].sum()

    # Participant support
    iso2pat = snvs[["iso_id","patient_id"]].drop_duplicates()
    support = (snv_feat[["iso_id","feature"]]
               .merge(iso2pat, on="iso_id", how="left")
               .groupby("feature", observed=True)["patient_id"].nunique())

    # Per-isolate totals for analytical null
    iso_totals = (gene_mass.groupby("iso_id", observed=True)["m_gene"]
                         .agg(T="sum", S2=lambda x: np.square(x).sum())
                         .astype(float))

    return obs, support, iso_totals, snv_feat, events

def iso_pq(fm: pd.DataFrame) -> tuple[dict, dict]:
    """
    Build per-isolate feature probabilities for the analytical null.
    For a random gene g ~ Uniform(genes-with-features in iso i),
    mass assigned to feature k is w_gk = 1/deg_g if k ∈ features(g) else 0.
      p_i,k = E[w_gk]  = (1/G_i) * Σ_g w_gk
      q_i,k = E[w_gk^2]= (1/G_i) * Σ_g w_gk^2
    """
    # how many genes-with-features per isolate (universe size)
    Gi = (fm.drop_duplicates(["iso_id","locus_tag"])
            .groupby("iso_id").size().astype(float).rename("Gi")).reset_index()

    # sum of weights and squared weights per (iso, feature)
    agg = (fm.assign(w2 = fm["w_split"]**2)
             .groupby(["iso_id","feature"], observed=True)
             .agg(sum_w=("w_split","sum"), sum_w2=("w2","sum"))
             .reset_index())

    # normalize by G_i to get unconditional expectations
    agg = agg.merge(Gi, on="iso_id", how="left")
    agg["p"] = agg["sum_w"] / agg["Gi"]
    agg["q"] = agg["sum_w2"] / agg["Gi"]

    # return as maps
    p_map = {iso: sub.set_index("feature")["p"] for iso, sub in agg.groupby("iso_id")}
    q_map = {iso: sub.set_index("feature")["q"] for iso, sub in agg.groupby("iso_id")}
    return p_map, q_map

def feature_enrichment_analytic(obs: pd.Series,
                                support: pd.Series,
                                iso_totals: pd.DataFrame,
                                p_map: dict, q_map: dict) -> pd.DataFrame:
    features = sorted(set(obs.index) |
                      set(f for d in p_map.values() for f in d.index))
    mu  = pd.Series(0.0, index=features)
    var = pd.Series(0.0, index=features)

    for iso, row in iso_totals.iterrows():
        Ti, S2i = row["T"], row["S2"]
        p = p_map.get(iso, pd.Series(0.0))
        q = q_map.get(iso, pd.Series(0.0))
        mu  = mu.add(Ti * p,                         fill_value=0.0)
        var = var.add(S2i * (q - p.pow(2)),          fill_value=0.0)

    var = var.clip(lower=1e-12)
    obs = obs.reindex(features, fill_value=0.0)
    enr = (obs + 1e-12) / (mu + 1e-12)
    z   = (obs - mu) / np.sqrt(var)
    p1  = 1.0 - norm.cdf(z)          # one-sided (excess mass)
    fdr = multipletests(p1, method="fdr_bh")[1]

    out = pd.DataFrame({
        "feature": features,
        "observed_mass": obs.values,
        "expected_mass": mu.values,
        "enrichment": enr.values,
        "z": z.values,
        "p_upper": p1,
        "FDR_upper": fdr,
        "n_participants": support.reindex(features, fill_value=0).values
    }).sort_values(["FDR_upper","p_upper"])
    return out

# Optional: tiny permutation fallback for very low-variance features
def permute_subset_pvalues(target_feats: list[str], B: int,
                           snv_feat: pd.DataFrame,      # from step 4 join (has iso_id, feature, mass)
                           fm: pd.DataFrame):           # feature map (iso_id, locus_tag, feature, w_split)
    import numpy.random as npr
    # build per-iso feature probabilities p_i,k from fm
    p_map, _ = iso_pq(fm)
    obs = snv_feat.groupby("feature")["mass"].sum()
    # per iso: list of SNV masses
    masses_by_iso = snv_feat.groupby("iso_id")["mass"].apply(list).to_dict()

    perm_exceed = {k:0 for k in target_feats}
    for _ in range(B):
        tot = {k:0.0 for k in target_feats}
        for iso, masses in masses_by_iso.items():
            p = p_map.get(iso, pd.Series(dtype=float))
            if p.empty: continue
            feats = p.index.to_numpy()
            probs = (p / p.sum()).to_numpy() if p.sum()>0 else None
            if probs is None: continue
            # draw a feature for each mass, add mass to that feature
            picks = feats[npr.choice(len(feats), size=len(masses), p=probs, replace=True)]
            for m, f in zip(masses, picks):
                if f in tot: tot[f] += m
        for k in target_feats:
            if tot.get(k,0.0) >= obs.get(k,0.0): perm_exceed[k]+=1

    return {k: (perm_exceed[k]+1)/(B+1) for k in target_feats}

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
    log.info("STARTING WEIGHTED SNV FEATURE ENRICHMENT ANALYSIS")
    log.info("="*60)
    log.info(f"Output directory: {args.out_dir}")
    log.info(f"Feature type: {args.feature_type}")
    log.info(f"Group filter: {args.group}")
    
    out_dir = args.out_dir; out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading input data...")
    meta = load_meta(args.meta_file)
    snv  = load_snvs(args.snv_path)
    mg   = load_mag_genes(args.mag_genes)
    flt  = load_filter(args.filter_file)
    
    # CRITICAL: Split locus tags containing '/' into separate rows and divide freq_range
    log.info("Splitting locus tags and dividing freq_range weights...")
    snv = split_locus_tags(snv, weight_col="freq_range")
    
    # Apply frequency floor filtering if specified
    if args.freq_floor > 0.0:
        log.info(f"Applying frequency floor filter (freq_range >= {args.freq_floor})...")
        before_count = len(snv)
        snv = snv[snv["freq_range"] >= args.freq_floor].copy()
        after_count = len(snv)
        log.info(f"Frequency floor filtering: {before_count} -> {after_count} SNVs ({after_count/before_count:.1%} retained)")

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

    # Build feature map from MAG genes
    log.info("Building feature mapping from MAG genes...")
    fm = build_feature_map(mg, args.feature_type)
    log.info(f"Built feature map with {len(fm)} gene-feature pairs")
    
    # Verify w_split weights are sane (no zero/NaN)
    assert (fm["w_split"] > 0).all(), "w_split must be > 0"

    # Performance optimization: use categoricals for groupbys/join keys
    log.info("Converting to categorical types for performance...")
    snv["iso_id"] = snv["iso_id"].astype("category")
    snv["locus_tag"] = snv["locus_tag"].astype("category")
    fm["iso_id"] = fm["iso_id"].astype("category")
    fm["locus_tag"] = fm["locus_tag"].astype("category")
    fm["feature"] = fm["feature"].astype("category")

    # Trim columns before merges to reduce copy cost
    snv = snv[["patient_id","bin","iso_id","locus_tag","freq_range"]].copy()

    # Coverage diagnostics: how many SNVs/mass map to any feature?
    log.info("Checking SNV-to-feature mapping coverage...")
    fm_pairs = fm[["iso_id","locus_tag"]].drop_duplicates()
    snv0 = snv.assign(m_base=snv["freq_range"].astype(float))
    snvJ = snv0.merge(fm_pairs, on=["iso_id","locus_tag"], how="inner", sort=False)

    rows_cov = len(snvJ) / max(1, len(snv0))
    mass_cov = snvJ["m_base"].sum() / max(1e-12, snv0["m_base"].sum())
    log.info(f"[COVERAGE] SNV rows with any feature: {rows_cov:.2%}; mass coverage: {mass_cov:.2%}")

    cov_iso = (snvJ.groupby("iso_id")["m_base"].sum() /
               snv0.groupby("iso_id")["m_base"].sum()).fillna(0.0)
    low = cov_iso[cov_iso < 0.5]
    if len(low):
        log.warning("[COVERAGE] Isolates with <50%% SNV mass mapped to features (n=%d): %s",
                    len(low), list(low.index[:10]))

    # Compute observed feature masses (all SNVs, weighted by freq_range)
    log.info("Computing observed feature masses from all SNVs...")
    obs, support, iso_totals, snv_feat, n_events = observed_feature_mass(
        snv, fm, collapse=args.collapse, freq_power=args.freq_power
    )
    iso_totals.index.name = "iso_id"  # (for the loop in feature_enrichment_analytic)
    log.info(f"Observed mass: {obs.sum():.3f} over {len(obs[obs > 0])} non-zero features")
    log.info(f"Top 10 observed features: {obs.nlargest(10).to_dict()}")

    # Compute null expectation analytically
    log.info("Computing analytical null expectation...")
    p_map, q_map = iso_pq(fm)
    
    # Post-condition checks (should all hold now)
    for iso, p in p_map.items():
        s = float(p.sum())
        if not np.isfinite(s) or abs(s - 1.0) > 1e-8:
            log.warning("[P-MAP] Sum_k p_{%s,k} = %.6f (should be 1).", iso, s)
    
    res = feature_enrichment_analytic(obs, support, iso_totals, p_map, q_map)

    # Mass conservation sanity checks
    sum_T  = float(iso_totals["T"].sum())
    sum_mu = float(res["expected_mass"].sum())
    sum_obs = float(res["observed_mass"].sum())
    log.info(f"[CHECK] ΣT={sum_T:.6f}  Σexpected={sum_mu:.6f}  Σobserved={sum_obs:.6f}")

    if abs(sum_mu - sum_T) > 1e-6 * max(1.0, sum_T):
        log.warning("[CHECK] Expected mass does NOT sum to ΣT. Inspect p_i,k construction.")
    if abs(sum_obs - sum_T) > 1e-6 * max(1.0, sum_T):
        log.warning("[CHECK] Observed mass does NOT equal ΣT. Inspect join/split logic.")

    # Add event counts to results and filter by minimum events
    res["n_events"] = n_events.reindex(res["feature"]).fillna(0).astype(int).values
    
    # Log event distribution
    log.info("[EVENTS] %d / %d features have ≥%d events",
             int((res["n_events"] >= args.min_events).sum()), len(res), args.min_events)
    
    # Apply the event-count filter before significance filtering
    res = res[res["n_events"] >= args.min_events].copy()
    log.info(f"After event filtering: {len(res)} features remain")

    # Filter significant results
    sig = res[(res["FDR_upper"] < args.p_thresh) &
              (res["enrichment"] >= args.min_enrichment) &
              (res["n_participants"] >= args.min_participants)].copy()

    stem = f"all_snv_weighted_enrichment_{args.feature_type}_{args.group.lower()}"
    out_file = out_dir / f"{stem}.csv"
    sig_file = out_dir / f"{stem}_significant.csv"
    
    res.to_csv(out_file, index=False)
    sig.to_csv(sig_file, index=False)
    
    log.info("Wrote: %s (%d features); significant: %s (%d)",
             out_file, len(res), sig_file, len(sig))

    with open(out_dir / "run_params.json","w") as fh:
        json.dump({
            "feature_type": args.feature_type,
            "group": args.group,
            "method": "analytical_weighted_snv",
            "p_thresh": args.p_thresh,
            "min_enrichment": args.min_enrichment,
            "min_participants": args.min_participants
        }, fh, indent=2)

def parse_args():
    ap = argparse.ArgumentParser("Weighted SNV feature enrichment (analytical null expectation)")
    ap.add_argument("--snv-path", type=Path, required=True,
                    help="Parquet/CSV file or directory with SNVs (patient_id, bin, freq_range, locus_tag)")
    ap.add_argument("--meta-file", type=Path, required=True,
                    help="Metadata with patient_id, diagnosis")
    ap.add_argument("--mag-genes", type=Path, required=True,
                    help="Compiled MAG gene TSV (patient_id, bin, locus_tag, cluster_rep, rep_gene, rep_dbxrefs)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--group", choices=["ALL","IBD","nonIBD"], default="IBD")
    ap.add_argument("--feature-type", choices=["kegg","go","ec","gene"], default="kegg")
    ap.add_argument("--filter-file", type=Path, help="Optional TSV with patient_id and bin to keep")
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--min-enrichment", type=float, default=2.0)
    ap.add_argument("--min-participants", type=int, default=3)
    ap.add_argument("--min-events", type=int, default=5,
                    help="Minimum # of iso×gene mutation events per feature to test")
    ap.add_argument("--freq-power", type=float, default=2.0,
                    help="Exponent α applied to per-site freq_range before gene collapse")
    ap.add_argument("--collapse", choices=["prob_any","max","sum"], default="prob_any",
                    help="How to collapse multiple SNVs per gene")
    ap.add_argument("--freq-floor", type=float, default=0.0,
                    help="Minimum freq_range to include (pre-filter low-frequency sites)")
    return ap.parse_args()

if __name__ == "__main__":
    run(parse_args())
