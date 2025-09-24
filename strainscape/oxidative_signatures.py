#!/usr/bin/env python3
"""
oxidative_signature_duckdb.py
Call exactly like the previous script – the API is unchanged.
Requires:  pip install duckdb pandas numpy statsmodels pyarrow
"""

import argparse
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -------------------------------------------------------------------- #
# 0. helper ­— vectorised BH FDR                                        #
# -------------------------------------------------------------------- #
def bh_adjust(pvals):
    return multipletests(pvals, method="fdr_bh")[1]


# -------------------------------------------------------------------- #
# 1.  read the small tables                                             #
# -------------------------------------------------------------------- #
def read_meta(path, group_col):
    cols = ["Participant ID", group_col]
    df = pd.read_csv(path, sep=None, engine="python").loc[:, cols]
    df = df.rename(columns={"Participant ID": "patient_id", group_col: "group"})
    if "group" not in df.columns:
        raise ValueError(f"Column '{group_col}' not present in metadata.")
    return df

def read_bin_identify(path):
    keep = ["user_genome", "classification", "patient_id"]
    return pd.read_csv(path, sep="\t").loc[:, keep]


# -------------------------------------------------------------------- #
# 2.  core work – all done inside DuckDB                                #
# -------------------------------------------------------------------- #
DUCK_SQL_TEMPLATE = """
WITH snv AS (
  SELECT
      patient_id,
      bin,
      slope,
      ref_base,
      new_base,
      freq_range,
      /* compute substitution on the fly */
      CASE WHEN slope < 0
           THEN new_base || '>' || ref_base
           ELSE ref_base || '>' || new_base
      END AS substitution
  FROM read_parquet('{snv_glob}')
  WHERE freq_range >= 0.7
),

joined AS (
  SELECT
      m."group",
      b.classification,
      s.patient_id,
      s.bin,
      s.substitution
  FROM snv s
  JOIN meta_df  m ON s.patient_id = m.patient_id
  JOIN bins_df  b ON s.bin       = b.user_genome
                 AND s.patient_id = b.patient_id
)

/* ------------------------------------------------------------------ */
/* sweep_counts – already aggregated; orders of magnitude smaller     */
SELECT
    "group",
    regexp_replace(classification, '.*s__', '') AS species,
    patient_id,
    bin,
    substitution,
    COUNT(*) AS count
FROM joined
GROUP BY 1,2,3,4,5
"""

def build_sweep_counts(con, snv_glob):
    query = DUCK_SQL_TEMPLATE.format(snv_glob=snv_glob)
    return con.execute(query).df()


# -------------------------------------------------------------------- #
# 3.  downstream pandas (unchanged but on far smaller tables)          #
# -------------------------------------------------------------------- #
def make_sweeps_prop(sweep_counts):
    all_subs = np.sort(sweep_counts["substitution"].unique())

    def _complete(g):
        g = g.set_index("substitution")
        g = g.reindex(all_subs, fill_value=0)
        g = g.reset_index()
        g["prop"] = g["count"] / g["count"].sum()
        return g

    sweeps_prop = (sweep_counts
        .groupby(["group", "species", "patient_id", "bin"], observed=True, as_index=False)
        .apply(_complete)
        .reset_index(drop=True)
    )
    return sweeps_prop


def pairwise_welch(df, grp):
    mask = df["group"].isin(["nonIBD", grp])
    keep_species = (df[mask]
                      .groupby(["species", "group"])
                      .size()
                      .unstack(fill_value=0))
    if grp in keep_species.columns and "nonIBD" in keep_species.columns:
        keep_species = keep_species[(keep_species["nonIBD"] >= 2) & (keep_species[grp] >= 2)].index
    else:
        keep_species = []

    rows = []
    for sp in keep_species:
        sub = df[(df["species"] == sp) & (df["group"].isin(["nonIBD", grp]))]

        nonibd = sub.loc[sub["group"] == "nonIBD", "prop_ox"]
        target = sub.loc[sub["group"] == grp,         "prop_ox"]

        stat, p = ttest_ind(target, nonibd, equal_var=False, nan_policy="omit")
        rows.append({"species": sp,
                     "comparison": f"{grp} vs nonIBD",
                     "estimate": target.mean() - nonibd.mean(),
                     "p_value": p})

    res = pd.DataFrame(rows)
    if not res.empty:
        res["p_adj"] = bh_adjust(res["p_value"].values)
    return res

# --- new helpers --------------------------------------------------------------
def make_global_prop(sweep_counts):
    """
    Collapse across species → counts per (patient_id, bin, group, substitution),
    then convert to proportions.
    """
    # 1. sum over species
    collapsed = (sweep_counts
        .groupby(["group", "patient_id", "bin", "substitution"], observed=True, as_index=False)
        .agg(count=("count", "sum"))
    )

    # 2. proportions within each patient/bin
    def _add_prop(g):
        g = g.copy()
        g["prop"] = g["count"] / g["count"].sum()
        return g

    global_prop = (collapsed
        .groupby(["group", "patient_id", "bin"], as_index=False, observed=True)
        .apply(_add_prop)
        .reset_index(drop=True)
    )
    return global_prop


def pairwise_welch_global(df):
    """Welch t-tests for every substitution (collapsed over species)."""
    rows = []
    for sub in df["substitution"].unique():
        sub_df = df[df["substitution"] == sub]

        for grp in [g for g in sub_df["group"].unique() if g != "nonIBD"]:
            nonibd = sub_df.loc[sub_df["group"] == "nonIBD", "prop"]
            target = sub_df.loc[sub_df["group"] == grp,      "prop"]

            # need at least 2 samples per arm
            if len(nonibd) >= 2 and len(target) >= 2:
                stat, p = ttest_ind(target, nonibd, equal_var=False, nan_policy="omit")
                rows.append({"substitution": sub,
                             "comparison": f"{grp} vs nonIBD",
                             "estimate": target.mean() - nonibd.mean(),
                             "p_value": p})

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj"] = bh_adjust(out["p_value"].values)
    return out


def run_once(meta, bins, snv_glob, out_prefix, output_dir):
    # ---------------- DuckDB phase ----------------------------------- #
    con = duckdb.connect(database=":memory:")
    con.register("meta_df", meta)
    con.register("bins_df", bins)

    sweep_counts = build_sweep_counts(con, snv_glob)
    con.close()

    # ---------------- pandas phase ----------------------------------- #
    sweeps_prop = make_sweeps_prop(sweep_counts)

    ox_sig = (sweeps_prop
        .query("substitution in ['C>A', 'G>T']")
        .groupby(["patient_id", "bin", "group", "species"], observed=True, as_index=False)
        .agg(prop_ox=("prop", "sum"))
    )

    comps   = [g for g in ox_sig["group"].unique() if g != "nonIBD"]
    stats   = pd.concat([pairwise_welch(ox_sig, g) for g in comps],
                        ignore_index=True) if comps else pd.DataFrame()

    # Ensure consistent data types for parquet writing
    sweeps_prop = sweeps_prop.astype({
        'group': 'string',
        'species': 'string', 
        'patient_id': 'string',
        'bin': 'string',
        'substitution': 'string',
        'count': 'int64',
        'prop': 'float64'
    })
    sweeps_prop.to_parquet(f"{output_dir}/{out_prefix}_SPECIES_transversion_prop.parquet", index=False)
    stats.to_csv       (f"{output_dir}/{out_prefix}_SPECIES_transversion_stats.tsv", sep="\t", index=False)

    global_prop = make_global_prop(sweep_counts)
    global_stats = pairwise_welch_global(global_prop)

    global_prop.to_parquet(f"{output_dir}/{out_prefix}_GLOBAL_transversion_prop.parquet", index=False)
    global_stats.to_csv(f"{output_dir}/{out_prefix}_GLOBAL_transversion_stats.tsv", sep="\t", index=False)


# -------------------------------------------------------------------- #
# 4.  CLI                                                              #
# -------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--snv_dir",      required=True)          # folder with *.parquet
    p.add_argument("--meta",         required=True)
    p.add_argument("--bin_identify", required=True)
    p.add_argument("--group_cols",   default="group",
                   help="Comma-separated columns in meta to test.")
    args = p.parse_args()

    bins = read_bin_identify(args.bin_identify)
    snv_glob = str(Path(args.snv_dir) / "*.parquet")         # DuckDB accepts globs

    output_dir = "/Users/reneeoles/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Research/Strain_Evolution/iHMP/output/transversions"
    for gcol in args.group_cols.split(","):
        meta = read_meta(args.meta, gcol)
        prefix = Path(args.snv_dir).stem + f"_{gcol}"
        run_once(meta, bins, snv_glob, prefix, output_dir)


if __name__ == "__main__":
    main()
