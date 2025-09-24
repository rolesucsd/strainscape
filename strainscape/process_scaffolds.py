#!/usr/bin/env python3
"""
ETL for genome/scaffold metrics (R→Python translation).

Inputs (TSV):
  --genome_info   : combined_genome_info.tsv  (inStrain genome table per sample)
  --checkm2       : checkm_merged.tsv         (must contain 'Bin Id', 'Completeness', 'Contamination')
  --taxonomy      : gtdbtk.bac120.summary.tsv (must contain 'user_genome', 'classification')
  --scaffolds     : combined_scaffold_info.tsv (inStrain scaffold table)

  Outputs (TSV in --outdir):
  - master_genomes_per_sample.tsv  (ALL genomes per sample with ALL columns and derived metrics & taxonomy)
  - master_genomes_per_sample_filtered.tsv  (QC-passed genomes per sample with derived metrics & taxonomy)
  - bin_aggregates_across_samples.tsv      (per-bin aggregate across samples)
  - genome_summary_across_samples.tsv      (per-sample aggregate of genome metrics and diversity)
  - scaffold_coverage_cv_by_bin_sample.tsv (per bin×sample: n_scaffolds, coverage CV, high flag)
  - scaffold_stats_detailed.tsv            (scaffold rows with bin_id/scaf_id/scaf_cov parsed)
"""

from __future__ import annotations
import argparse, sys, re, os, math, logging
from pathlib import Path
from typing import Optional, Sequence, Dict

import numpy as np
import pandas as pd

# ----------------------------- #
#            Loggers            #
# ----------------------------- #
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

L = logging.getLogger("etl")

# ----------------------------- #
#       Robust TSV reader       #
# ----------------------------- #
def read_tsv_robust(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

def to_float(s: pd.Series, na=np.nan) -> pd.Series:
    try:
        return pd.to_numeric(s, errors="coerce")
    except Exception:
        return pd.Series([na]*len(s))

def to_int(s: pd.Series, na=np.nan) -> pd.Series:
    try:
        v = pd.to_numeric(s, errors="coerce")
        # keep NaN for non-integers
        return v.astype("Int64")
    except Exception:
        return pd.Series([na]*len(s)).astype("Int64")

def has(df: Optional[pd.DataFrame], col: str) -> bool:
    return (df is not None) and (col in df.columns)

# ----------------------------- #
#     GTDB classification →     #
#        taxonomic ranks        #
# ----------------------------- #
RANKS = ["domain","phylum","class","order","family","genus","species"]
RANK_PREFIX = re.compile(r"(?i)(?:^|;)\s*[dkpcofgs]__")

def strip_rank_tags(s: str) -> str:
    # remove d__/k__/D_0__ etc.
    # Replace the regex pattern to properly remove rank prefixes
    return re.sub(r'[dkpcofgs]__', '', s)

def parse_gtdb_taxonomy(tax_df: pd.DataFrame, col: str="classification") -> pd.DataFrame:
    if tax_df is None or col not in tax_df.columns:
        return tax_df if tax_df is not None else pd.DataFrame()
    tx = tax_df.copy()
    tx[".tax"] = tx[col].astype(str).map(strip_rank_tags)
    parts = tx[".tax"].str.split(";", expand=True)
    parts = parts.iloc[:, :len(RANKS)].rename(columns=dict(enumerate(RANKS)))
    parts = parts.apply(lambda c: c.str.strip().replace({"": np.nan}), axis=0)
    tx = pd.concat([tx.drop(columns=[".tax"], errors="ignore"), parts], axis=1)
    return tx

# ----------------------------- #
#      Harmonize genome_info    #
# ----------------------------- #
def harmonize_genome_info(gi: pd.DataFrame) -> pd.DataFrame:
    df = gi.copy()

    # Ensure Sample exists as str
    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str)
    else:
        df["Sample"] = np.nan

    # Column aliases
    colmap = {
        "coverage_median": "median_coverage",
        "length": "genome_length",
        "SNV_count": "SNVs",
        "SNS_count": "SNSs",
        "nucl_diversity": "nucleotide_diversity",
        "consensus_divergent_sites": "consensus_div_sites",
        "population_divergent_sites": "pop_div_sites",
        "genome": "bin_id",
        "breadth": "breadth",  # already named ok in most exports
    }

    # Create target columns if present in source
    for src, dst in colmap.items():
        if src in df.columns:
            df[dst] = df[src]

    # Coerce numerics safely
    num_cols = [
        "median_coverage","genome_length","SNVs","SNSs","nucleotide_diversity",
        "consensus_div_sites","pop_div_sites","breadth"
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = to_float(df[c])

    # Keep bin_id even if named differently
    if "bin_id" not in df.columns:
        # Try alternatives
        for alt in ("user_genome","genome","bin","Bin Id","bin"):
            if alt in df.columns:
                df["bin_id"] = df[alt].astype(str)
                break

    # Standardize columns presence
    for req in ["bin_id","Sample","genome_length","median_coverage","breadth"]:
        if req not in df.columns:
            df[req] = np.nan

    return df

# ----------------------------- #
#         CheckM2 & Tax         #
# ----------------------------- #
def harmonize_checkm2(chk: Optional[pd.DataFrame]) -> pd.DataFrame:
    if chk is None or len(chk)==0:
        return pd.DataFrame(columns=["bin_id","completeness","contamination"])
    df = chk.copy()
    # Expect 'Bin Id', 'Completeness', 'Contamination'
    if "Bin Id" in df.columns:
        df["bin_id"] = df["Bin Id"].astype(str)
    elif "bin_id" not in df.columns:
        # last resort: try 'bin' or 'genome'
        for alt in ("bin","genome","user_genome"):
            if alt in df.columns:
                df["bin_id"] = df[alt].astype(str)
                break
    df["completeness"] = to_float(df.get("Completeness", pd.Series([np.nan]*len(df))))
    df["contamination"] = to_float(df.get("Contamination", pd.Series([np.nan]*len(df))))
    return df[["bin_id","completeness","contamination"]].drop_duplicates()

def harmonize_taxonomy(tax: Optional[pd.DataFrame]) -> pd.DataFrame:
    if tax is None or len(tax)==0:
        return pd.DataFrame(columns=["bin_id"] + RANKS)
    tx = tax.copy()
    if "user_genome" in tx.columns:
        tx["bin_id"] = tx["user_genome"].astype(str)
    elif "bin_id" not in tx.columns:
        # attempt fallbacks
        for alt in ("bin","genome","Bin Id"):
            if alt in tx.columns:
                tx["bin_id"] = tx[alt].astype(str)
                break
    tx = parse_gtdb_taxonomy(tx, "classification")
    keep = ["bin_id"] + [r for r in RANKS if r in tx.columns]
    return tx[keep].drop_duplicates()

# ----------------------------- #
#  Scaffold parsing + coverage  #
# ----------------------------- #
def parse_scaffolds(scf: Optional[pd.DataFrame], cov_threshold: float = 10.0, breadth_threshold: float = 0.9) -> pd.DataFrame:
    """
    Expect columns: scaffold, Sample, and either coverage_median or coverage, and length/breadth optionally.
    Attempt to derive bin_id, scaf_id from 'scaffold' if not provided:
      e.g., 'binX.fa|NODE_12' -> bin_id='binX', scaf_id='NODE_12'
    
    Filters out individual scaffolds that don't meet coverage and breadth thresholds.
    """
    if scf is None or len(scf)==0:
        return pd.DataFrame(columns=["Sample","bin_id","scaf_id","scaffold","scaf_cov","length","breadth"])

    df = scf.copy()

    for c in ("Sample","scaffold"):
        if c not in df.columns:
            df[c] = np.nan

    # Coverage selection
    if "coverage_median" in df.columns:
        df["scaf_cov"] = to_float(df["coverage_median"])
    elif "coverage" in df.columns:
        df["scaf_cov"] = to_float(df["coverage"])
    else:
        df["scaf_cov"] = np.nan

    # length / breadth optional
    if "length" in df.columns:
        df["length"] = to_float(df["length"])
    else:
        df["length"] = np.nan
    if "breadth" in df.columns:
        df["breadth"] = to_float(df["breadth"])
    else:
        df["breadth"] = np.nan

    # bin_id / scaf_id parsing
    if "bin_id" not in df.columns:
        # Try to parse from 'scaffold'
        # pattern: "<bin>(.fa|.fasta)?|<scaf>"
        bin_id = []
        scaf_id = []
        for s in df["scaffold"].astype(str):
            if "|" in s:
                left, right = s.split("|", 1)
                left = re.sub(r"\.(fa|fasta|fna|fas)$", "", left)
                bin_id.append(left)
                scaf_id.append(right)
            else:
                bin_id.append(np.nan)
                scaf_id.append(s)
        df["bin_id"] = bin_id
        df["scaf_id"] = scaf_id
    else:
        if "scaf_id" not in df.columns:
            df["scaf_id"] = df["scaffold"].astype(str)

    # Clean types
    df["Sample"] = df["Sample"].astype(str)
    df["bin_id"]  = df["bin_id"].astype(str)

    # Filter individual scaffolds based on coverage and breadth thresholds
    initial_count = len(df)
    
    # Filter by coverage threshold
    if "scaf_cov" in df.columns:
        df = df[df["scaf_cov"] >= cov_threshold]
        coverage_filtered = initial_count - len(df)
        if coverage_filtered > 0:
            L.info(f"Scaffold coverage filter: removed {coverage_filtered:,} scaffolds below {cov_threshold}x coverage")
    
    # Filter by breadth threshold (if available)
    if "breadth" in df.columns and not df["breadth"].isna().all():
        before_breadth = len(df)
        df = df[df["breadth"] >= breadth_threshold]
        breadth_filtered = before_breadth - len(df)
        if breadth_filtered > 0:
            L.info(f"Scaffold breadth filter: removed {breadth_filtered:,} scaffolds below {breadth_threshold} breadth")
    
    final_count = len(df)
    if final_count < initial_count:
        L.info(f"Individual scaffold filtering: {initial_count:,} → {final_count:,} scaffolds (removed {initial_count - final_count:,})")

    cols = ["Sample","bin_id","scaf_id","scaffold","scaf_cov","length","breadth"]
    return df[cols]

def scaffold_cv_by_bin_sample(scf_parsed: pd.DataFrame, cv_threshold: float) -> pd.DataFrame:
    """
    Returns per bin×Sample:
      n_scaffolds, scaf_cov_cv, high_scaffold_cov_cv (bool)
    """
    if scf_parsed is None or len(scf_parsed)==0:
        return pd.DataFrame(columns=["Sample","bin_id","n_scaffolds","scaf_cov_cv","high_scaffold_cov_cv"])
    df = scf_parsed.dropna(subset=["bin_id","Sample"]).copy()
    grp = df.groupby(["bin_id","Sample"], dropna=False)
    agg = grp["scaf_cov"].agg(["count","mean","std"]).reset_index()
    agg.rename(columns={"count":"n_scaffolds","std":"sd","mean":"mean"}, inplace=True)
    # CV = sd / mean (guard against zero/NaN)
    agg["scaf_cov_cv"] = np.where((agg["mean"].abs() > 0) & np.isfinite(agg["sd"]),
                                  agg["sd"] / agg["mean"].abs(),
                                  np.nan)
    agg["high_scaffold_cov_cv"] = (agg["scaf_cov_cv"] >= float(cv_threshold))
    # Where CV is NaN (e.g., n=1), mark flag False
    agg.loc[~np.isfinite(agg["scaf_cov_cv"]), "high_scaffold_cov_cv"] = False
    return agg[["Sample","bin_id","n_scaffolds","scaf_cov_cv","high_scaffold_cov_cv"]]

# ----------------------------- #
#       Core ETL pipeline       #
# ----------------------------- #
def run_etl(
    genome_info: str,
    checkm2: str,
    taxonomy: str,
    scaffolds: str,
    outdir: str,
    cov_threshold: float,
    breadth_threshold: float,
    comp_threshold: float,
    cont_threshold: float,
    cv_scaf_threshold: float,
):
    Path(outdir).mkdir(parents=True, exist_ok=True)

    gi_raw = read_tsv_robust(genome_info)
    chk_raw = read_tsv_robust(checkm2)
    tax_raw = read_tsv_robust(taxonomy)
    scf_raw = read_tsv_robust(scaffolds)

    if gi_raw is None or len(gi_raw)==0:
        L.error("genome_info is empty or missing.")
        sys.exit(2)

    gi = harmonize_genome_info(gi_raw)
    chk = harmonize_checkm2(chk_raw)
    tax = harmonize_taxonomy(tax_raw)
    scf = parse_scaffolds(scf_raw, cov_threshold=cov_threshold, breadth_threshold=breadth_threshold)

    L.info(f"Loaded: genome_info={len(gi):,}, checkm2={len(chk):,}, taxonomy={len(tax):,}, scaffolds={len(scf):,}")

    # Join CheckM2 + Taxonomy
    df = gi.merge(chk, on="bin_id", how="left")
    df = df.merge(tax, on="bin_id", how="left")

    # Prevalence calculation across samples
    medcov = to_float(df.get("median_coverage", pd.Series(np.nan, index=df.index)))
    breadth = to_float(df.get("breadth", pd.Series(np.nan, index=df.index)))
    df["is_present"] = (np.isfinite(medcov) & (medcov >= cov_threshold) &
                        np.isfinite(breadth) & (breadth >= breadth_threshold))

    prevalence = (
        df.groupby("bin_id", dropna=False)["is_present"]
          .sum(min_count=1)
          .fillna(0)
          .astype(int)
          .rename("prevalence")
          .reset_index()
    )
    df = df.merge(prevalence, on="bin_id", how="left")
    df["prevalence"] = df["prevalence"].fillna(0).astype(int)

    # Coverage-aware denominators + per-Mbp rates
    df["covered_bp"] = np.where(
        np.isfinite(df["breadth"]) & np.isfinite(df["genome_length"]),
        df["breadth"] * df["genome_length"],
        np.nan
    )

    def per_mbp(numer, denom):
        ok = np.isfinite(numer) & np.isfinite(denom) & (denom > 0)
        out = np.full(len(df), np.nan, dtype=float)
        out[ok] = (numer[ok] / denom[ok]) * 1e6
        return out

    for col in ("SNVs","SNSs","consensus_div_sites","pop_div_sites"):
        if col not in df.columns:
            df[col] = np.nan

    df["snvs_per_mbp"] = per_mbp(to_float(df["SNVs"]), to_float(df["genome_length"]))
    df["snss_per_mbp"] = per_mbp(to_float(df["SNSs"]), to_float(df["genome_length"]))

    df["consensus_div_per_mbp"] = per_mbp(to_float(df["consensus_div_sites"]), to_float(df["genome_length"]))
    df["pop_div_per_mbp"]       = per_mbp(to_float(df["pop_div_sites"]),        to_float(df["genome_length"]))

    def per_covbp(numer, covered):
        ok = np.isfinite(numer) & np.isfinite(covered) & (covered > 0)
        out = np.full(len(df), np.nan, dtype=float)
        out[ok] = numer[ok] / covered[ok]
        return out

    df["consensus_div_per_covered_bp"] = per_covbp(to_float(df["consensus_div_sites"]), to_float(df["covered_bp"]))
    df["pop_div_per_covered_bp"]       = per_covbp(to_float(df["pop_div_sites"]),        to_float(df["covered_bp"]))

    # QC gates
    df["is_high_cov"] = (np.isfinite(medcov) & (medcov >= cov_threshold) &
                         np.isfinite(breadth) & (breadth >= breadth_threshold))
    df["is_high_quality"] = (np.isfinite(df.get("completeness")) & (to_float(df["completeness"]) >= comp_threshold) &
                             np.isfinite(df.get("contamination")) & (to_float(df["contamination"]) <= cont_threshold))

    comp = to_float(df.get("completeness", pd.Series(np.nan, index=df.index)))
    cont = to_float(df.get("contamination", pd.Series(np.nan, index=df.index)))
    df["quality_score"] = np.where(np.isfinite(comp) & np.isfinite(cont), comp - 5.0*cont, np.nan)

    # -------- Scaffold stats (CV) --------
    scf_cv = scaffold_cv_by_bin_sample(scf, cv_scaf_threshold)

    # Optionally merge scaffold CV features into main df (per bin×Sample)
    df = df.merge(scf_cv, on=["bin_id","Sample"], how="left")

    # ----------------------------- #
    #     Write master tables       #
    # ----------------------------- #
    keep_tax = [c for c in (
        "Domain","Phylum","Class","Order","Family","Genus","Species",
        "domain","phylum","class","order","family","genus","species"
    ) if c in df.columns]

    # Fix R naming typos: use *_per_covered_bp consistently
    select_cols = [
        "Sample","bin_id",
        "genome_length",
        "median_coverage","breadth","covered_bp",
        "SNVs","SNSs",
        "nucleotide_diversity",
        "consensus_div_sites","pop_div_sites",
        "snvs_per_mbp","snss_per_mbp",
        "consensus_div_per_mbp","pop_div_per_mbp",
        "consensus_div_per_covered_bp","pop_div_per_covered_bp",
        "prevalence",
        "completeness","contamination","quality_score",
        "is_high_cov","is_high_quality",
        "scaf_cov_cv","n_scaffolds","high_scaffold_cov_cv",
    ] + keep_tax

    # Write comprehensive master table with ALL columns and ALL rows
    df.to_csv(Path(outdir, "master_genomes_per_sample.tsv"), sep="\t", index=False)
    
    # Write filtered master table with only QC-passed genomes
    df_filter = df.loc[df["is_high_quality"] & df["is_high_cov"], select_cols].drop_duplicates()
    df_filter.to_csv(Path(outdir, "master_genomes_per_sample_filtered.tsv"), sep="\t", index=False)

    # Aggregate per bin across samples (only QC-passed)
    agg_src = df[df["is_high_quality"] & df["is_high_cov"]].copy()

    def first_nonnull(series: pd.Series):
        for v in series:
            if pd.notnull(v):
                return v
        return np.nan

    agg = (
        agg_src
        .groupby("bin_id", dropna=False)
        .agg(
            n_samples=("Sample", lambda s: s.dropna().nunique()),
            prevalence=("prevalence", "max"),
            med_cov=("median_coverage", "median"),
            med_breadth=("breadth", "median"),
            med_snvs_per_mbp=("snvs_per_mbp", "median"),
            med_pi=("nucleotide_diversity", "median"),
            genome_length=("genome_length", first_nonnull),
            phylum=("phylum", first_nonnull) if "phylum" in agg_src.columns else ("bin_id", "size"),
            genus=("genus", first_nonnull) if "genus" in agg_src.columns else ("bin_id", "size"),
            species=("species", first_nonnull) if "species" in agg_src.columns else ("bin_id", "size"),
        )
        .reset_index()
    )

    # If taxonomy missing, drop the placeholder tuple columns
    for c in ("phylum","genus","species"):
        if c in agg.columns and isinstance(agg[c].dtype, pd.CategoricalDtype):
            agg[c] = agg[c].astype(str)
        # Correct accidental multi-function result
        if agg[c].dtype == object and agg[c].apply(lambda x: isinstance(x, tuple)).any():
            agg.drop(columns=[c], inplace=True, errors="ignore")

    agg.to_csv(Path(outdir, "bin_aggregates_across_samples.tsv"), sep="\t", index=False)

    # Filter scaffolds to only include those from QC-passed bins
    qc_passed_bins = set(agg_src["bin_id"].unique())
    scf_filtered = scf[scf["bin_id"].isin(qc_passed_bins)].copy()
    scf_cv_filtered = scf_cv[scf_cv["bin_id"].isin(qc_passed_bins)].copy()
    
    L.info(f"Scaffold filtering: {len(scf):,} → {len(scf_filtered):,} scaffolds (kept only from QC-passed bins)")
    L.info(f"Scaffold CV filtering: {len(scf_cv):,} → {len(scf_cv_filtered):,} bin×sample combinations")

    # Write scaffold CV + detailed scaffold stats (filtered to QC-passed bins only)
    scf_cv_filtered.sort_values(["bin_id","Sample"]).to_csv(Path(outdir, "scaffold_coverage_cv_by_bin_sample.tsv"), sep="\t", index=False)
    scf_filtered.sort_values(["bin_id","Sample","scaf_id"]).to_csv(Path(outdir, "scaffold_stats_detailed.tsv"), sep="\t", index=False)

    L.info("Wrote outputs to %s", outdir)

# ----------------------------- #
#             CLI               #
# ----------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Genome/scaffold ETL (R→Python) with QC gates, prevalence, taxonomy, and scaffold CV."
    )
    p.add_argument("--genome_info",   default="input/combined_genome_info.tsv")
    p.add_argument("--checkm2",       default="input/checkm_merged.tsv")
    p.add_argument("--taxonomy",      default="input/gtdbtk.bac120.summary.tsv")
    p.add_argument("--scaffolds",     default="input/combined_scaffold_info.tsv")
    p.add_argument("--outdir",        default="output/tmp")
    p.add_argument("--cov_threshold",       type=float, default=10.0)
    p.add_argument("--breadth_threshold",   type=float, default=0.9)
    p.add_argument("--comp_threshold",      type=float, default=80.0)
    p.add_argument("--cont_threshold",      type=float, default=5.0)
    p.add_argument("--cv_scaf_threshold",   type=float, default=0.8)
    return p

def main(argv: Optional[Sequence[str]]=None):
    setup_logging()
    args = build_parser().parse_args(argv)
    run_etl(
        genome_info=args.genome_info,
        checkm2=args.checkm2,
        taxonomy=args.taxonomy,
        scaffolds=args.scaffolds,
        outdir=args.outdir,
        cov_threshold=args.cov_threshold,
        breadth_threshold=args.breadth_threshold,
        comp_threshold=args.comp_threshold,
        cont_threshold=args.cont_threshold,
        cv_scaf_threshold=args.cv_scaf_threshold,
    )

if __name__ == "__main__":
    main()
