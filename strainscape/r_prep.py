#!/usr/bin/env python3
import logging, math
from pathlib import Path
from typing import Optional, Tuple, Dict, Union, Sequence

import pandas as pd, numpy as np
import pyarrow.dataset as ds
import pyarrow.feather as feather
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(asctime)s %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Paths (edit as needed)
# ─────────────────────────────────────────────────────────────────────────────
DIR   = Path("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output")
META  = Path("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/metadata/hmp2_metadata_2018-08-20.csv")
BINS  = Path("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/patient_bin_summary.csv")
PREP  = Path("/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/prep")
PREP.mkdir(exist_ok=True)

# ───────────── helpers ─────────────

def clean_colnames(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = (df.columns.str.strip().str.lower()
                           .str.replace(r"[^0-9a-zA-Z]+", "_", regex=True)
                           .str.replace(r"_+", "_", regex=True)
                           .str.strip("_"))
    return df

def write_delta_histogram(
    snvs: pd.DataFrame,
    out_path: Path,
    edges: Sequence[float] = np.linspace(0, 1, 21)  # 0.05-wide bins by default
) -> None:
    """
    Compute |max_freq – min_freq| for each SNV, bin the values,
    and write a tidy table to `out_path` with columns:
        bin, count
    Ready for ggplot2::geom_col() or geom_bar().
    """
    # 1. ensure frequencies are within [0, 1]
    clipped = snvs[["min_freq", "max_freq"]].clip(lower=0, upper=1)

    # 2. delta = absolute difference
    delta = (clipped["max_freq"] - clipped["min_freq"]).abs()

    # 3. bin the deltas
    labels = [f"{edges[i]:.2f}–{edges[i+1]:.2f}" for i in range(len(edges) - 1)]
    binned = pd.cut(delta, bins=edges, labels=labels, right=False)

    # 4. histogram
    hist = (
        binned.value_counts(sort=False, dropna=False)
              .reset_index(name="count")
              .rename(columns={"index": "bin"})
              .sort_values("bin")
    )

    hist.to_csv(out_path, index=False)
    logger.info("Δ-frequency histogram written → %s", out_path)


def middle_slash_split(val: str) -> Optional[Tuple[str, str]]:
    parts = val.split("/")
    if len(parts) < 2:
        return None
    mid = math.ceil(len(parts) / 2)
    return "/".join(parts[:mid]), "/".join(parts[mid:])


def read_metadata(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(meta_path)
    meta = clean_colnames(meta)
    meta = (meta.rename(columns={"participant_id": "patient_id"})
                .assign(group=lambda d: d["diagnosis"].map({"nonIBD":"nonIBD","UC":"UC","CD":"CD"}))
                [["patient_id","group"]].drop_duplicates())
    logger.info("Metadata → %d patients", meta.patient_id.nunique())
    return meta


def summarise_bins(bins_path: Path, meta: pd.DataFrame) -> pd.DataFrame:
    bins = pd.read_csv(bins_path)
    num_cols = bins.select_dtypes(include=[np.number]).columns
    bins_sum = (bins.groupby(["patient_id","bin"], as_index=False)[num_cols].mean()
                     .merge(meta, on="patient_id", how="left"))
    bins_sum["group"] = pd.Categorical(bins_sum["group"],
                                        categories=["nonIBD","UC","CD"], ordered=True)
    logger.info("Bins summary: %d×%d", *bins_sum.shape)
    return bins_sum


def load_snvs_parquet(base_dir: Path) -> pd.DataFrame:
    table = ds.dataset(base_dir, format="parquet").to_table()
    snvs  = table.to_pandas(split_blocks=True, self_destruct=True)
    logger.info("SNVs loaded: %d×%d", *snvs.shape)
    return snvs


def compute_sweeps(snvs: pd.DataFrame) -> pd.DataFrame:
    snvs = snvs.copy()
    snvs["is_sweep"] = (snvs["min_freq"] - snvs["max_freq"]).abs() >= 0.7
    logger.info("Sweeps: %d (%.2f%%)", snvs.is_sweep.sum(), 100*snvs.is_sweep.mean())
    return snvs


def snvs_per_bin(snvs: pd.DataFrame) -> pd.DataFrame:
    counts = (snvs.groupby(["patient_id","bin","is_sweep"], as_index=False)
                   .size().pivot(index=["patient_id","bin"], columns="is_sweep", values="size")
                   .fillna(0).reset_index().rename(columns={True:"TRUE",False:"FALSE"}))
    counts.columns.name = None
    logger.info("SNV counts per bin: %d", len(counts))
    return counts


def filter_bins_summary(bins_sum: pd.DataFrame, snv_counts: pd.DataFrame) -> pd.DataFrame:
    df = bins_sum.merge(snv_counts, on=["patient_id","bin"], how="inner")
    before = len(df)
    feather.write_feather(df, PREP/"bins_summary_filter.feather")
    df = df[df["TRUE"] <= 500].dropna(subset=["TRUE","FALSE"])
    logger.info("Filter bins TRUE≤1000: %d → %d", before, len(df))
    df["prop"]  = df["TRUE"] / df["FALSE"].replace(0, np.nan)
    df["total"] = df["TRUE"] + df["FALSE"]
    return df

# ───────── enrichment helpers ─────────
_ArrayLike = Union[Sequence, np.ndarray, pd.Series]

def _is_all_nan(arr: _ArrayLike) -> bool:
    """True if every element is NA; works for lists/arrays."""
    return all(pd.isna(x) for x in arr)


def explode_feature(df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        val = row[feature_col]
        # Handle list/array cells from Arrow → pandas conversion
        if isinstance(val, (list, tuple, np.ndarray, pd.Series)):
            if len(val) == 0 or _is_all_nan(val):
                continue
            val = ",".join(str(x) for x in val if not pd.isna(x))
        elif pd.isna(val):
            continue

        val = str(val)
        mut_type = row.get("mutation_type")

        # intergenic middle‑slash
        if mut_type == "Intergenic" and "/" in val:
            halves = middle_slash_split(val)
            if halves:
                for half in halves:
                    new = row.copy(); new[feature_col] = half.strip(); rows.append(new)
                continue

        # comma split for "*term" columns
        if "term" in feature_col and "," in val:
            for part in val.split(','):
                new = row.copy(); new[feature_col] = part.strip(); rows.append(new)
            continue

        row[feature_col] = val.strip(); rows.append(row)

    expanded = pd.DataFrame(rows)
    expanded[feature_col] = expanded[feature_col].str.strip()
    expanded = expanded[~expanded[feature_col].isin(["", "nan", "hypothetical protein"])]
    return expanded


# ── replace the whole function ───────────────────────────────────────────────
def fisher_enrich_py(df: pd.DataFrame,
                     feature_col: str,
                     min_prop: float = 0.02,
                     min_abs:  int   = 0) -> pd.DataFrame:
    """
    One-sided Fisher enrichment
    """
    # explode nested fields, one row per SNV/feature
    expanded = explode_feature(df, feature_col)

    # drop rare features
    total_isolates = expanded[["patient_id", "bin"]].drop_duplicates().shape[0]
    feat_isolates  = (
        expanded.groupby(feature_col)
             .apply(lambda d: d[["patient_id", "bin"]]
                       .drop_duplicates()
                       .shape[0], include_groups=False)
             .reset_index(name="iso_n")
    )
    min_needed = max(min_abs, int(min_prop * total_isolates))
    keep_feats = feat_isolates.loc[feat_isolates.iso_n >= min_needed, feature_col]
    if keep_feats.empty:
        return pd.DataFrame()                       # nothing to test
    expanded = expanded[expanded[feature_col].isin(keep_feats)]

    group_sizes = (
        expanded[["patient_id", "bin", "group"]]
             .drop_duplicates()
             .groupby("group").size().to_dict()
    )
    feature_counts = expanded.groupby([feature_col, "group"], as_index=False).size()

    # Fisher for each comparison
    recs, pairs = [], {
        "CD_vs_nonIBD": ("CD",  "nonIBD"),
        "CD_vs_UC":     ("CD",  "UC"),
        "UC_vs_nonIBD": ("UC",  "nonIBD")
    }
    for cname, (g1, g2) in pairs.items():
        sub  = feature_counts[feature_counts.group.isin([g1, g2])]
        wide = sub.pivot(index=feature_col,
                         columns="group",
                         values="size").fillna(0)
        for feat, row in wide.iterrows():
            n1, n2 = int(row.get(g1, 0)), int(row.get(g2, 0))
            t1, t2 = int(group_sizes.get(g1, 0)), int(group_sizes.get(g2, 0))
            try:
                odds, p = fisher_exact([[n1, t1 - n1],
                                        [n2, t2 - n2]])
            except ValueError:                       # zero marginal
                odds, p = np.nan, 1.0               # neutral
            enr = ((n1 + 0.5) / (t1 + 1e-9)) / ((n2 + 0.5) / (t2 + 1e-9))
            recs.append(dict(comparison=cname, feature=feat,
                             **{f"n_{g1}": n1, f"total_{g1}": t1,
                                f"n_{g2}": n2, f"total_{g2}": t2},
                             enrichment=enr, odds_ratio=odds, p_value=p))
    out = pd.DataFrame.from_records(recs)
    if out.empty:
        return out

    # Benjamini–Hochberg (safe-fill NAs with 1.0)
    pvals = out.p_value.fillna(1.0).to_numpy()
    out["FDR"] = multipletests(pvals, method="fdr_bh")[1]

    # convenience columns for volcanoes
    out["logE"] = np.log10(out.enrichment.replace(0, np.nan))
    out["logP"] = -np.log10(out.p_value.replace(0, np.nan))
    out["sig"]  = (out.enrichment > 2) & (out.p_value < 0.05)
    return out


# ───────── main ─────────

def main():
    meta      = read_metadata(META)
    bins_sum  = summarise_bins(BINS, meta)
    feather.write_feather(bins_sum, PREP/"bins_summary.feather")

    logger.info("Loading SNVs parquet – grab coffee…")
    snvs_raw = load_snvs_parquet(DIR / "all_snvs")

    write_delta_histogram(
        snvs_raw,
        PREP / "freq_histogram.csv",
        edges=np.linspace(0, 1, 21)
    )

    snvs     = compute_sweeps(snvs_raw)
    sweeps    = snvs[snvs.is_sweep].copy()
    feather.write_feather(sweeps, PREP/"sweeps.feather")
    counts    = snvs_per_bin(snvs)
    feather.write_feather(counts, PREP/"snvs_per_bin.feather")

    bins_filt = filter_bins_summary(bins_sum, counts)

    sweeps_f  = sweeps.merge(bins_filt[["patient_id","bin"]])
    logger.info("Sweeps ≤1k: %d rows", len(sweeps_f))
    feather.write_feather(sweeps_f, PREP/"sweeps_filtered.feather")

    subsets = {
        "snv_sweeps": sweeps,
        "snv_sweeps_nosilent": sweeps[sweeps.mutation_type != "Silent"],
        "snv_sweeps_1k": sweeps_f,
        "snv_sweeps_1k_nosilent": sweeps_f[sweeps_f.mutation_type != "Silent"]
    }
    feat_cols = {"gene":"gene","product":"product","kegg_terms":"kegg_terms",
                 "go_terms":"go_terms","ec_terms":"ec_terms"}

    all_res = []
    for sname,df in subsets.items():
        df = df.merge(meta, on="patient_id", how="left")
        # Convert numpy arrays/lists to comma-separated strings before dropping duplicates
        def flatten_val(x):
            if isinstance(x, (list, np.ndarray, pd.Series)):
                # Remove NaNs and convert to string
                vals = [str(i) for i in x if not pd.isna(i)]
                return ','.join(vals)
            return x
        for col in ["gene", "product", "kegg_terms", "go_terms", "ec_terms"]:
            if col in df.columns:
                df[col] = df[col].apply(flatten_val)
        df = df.drop_duplicates(subset=["patient_id", "bin", "gene", "product", 
                                        "kegg_terms", "go_terms", "ec_terms"])
        logger.info("Subset %s (%d rows)", sname, len(df))
        for fname,col in feat_cols.items():
            if col not in df.columns: continue
            out = fisher_enrich_py(df, col)
            logger.info("  %s – %s: %d feats", sname, fname, len(out))
            if out.empty: continue
            out["FDR"] = multipletests(out.p_value, method="fdr_bh")[1]
            out["logE"] = np.log10(out.enrichment.replace(0, np.nan))
            out["logP"] = -np.log10(out.p_value.replace(0, np.nan))
            out["sig"]  = (out.enrichment>2) & (out.p_value<0.05)
            out.insert(0,"subset",sname); out.insert(1,"feature_type",fname)
            all_res.append(out)

    if not all_res:
        logger.warning("No enrichment results – abort")
        return

    res_df = pd.concat(all_res, ignore_index=True)
    feather.write_feather(res_df, PREP/"enrichment_results.feather")
    res_df.to_csv(PREP/"enrichment_results.csv", index=False)
    logger.info("All done – results in %s", PREP.resolve())

if __name__ == "__main__":
    main()
