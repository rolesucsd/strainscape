#!/usr/bin/env python3
"""
Process inStrain scaffold_info + STB + sample metadata.

Steps
─────
1.  Build one small DF mapping scaffold ↦ bin (from *.fa in bin_dir)
2.  Merge with STB  (adds Patient-/Sample-ID)
3.  Merge with scaffold_info.tsv, apply length / coverage / breadth filters
4.  Merge once with metadata (selected columns only)
5.  Write filtered table + log summary stats
"""

from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse, logging, sys, textwrap

# ────────────────── logging ──────────────────
def setup_logger(logf):
    fmt = "%(asctime)s %(levelname)s  %(message)s"
    h = logging.FileHandler(logf)
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[h, logging.StreamHandler(sys.stderr)])
    return logging.getLogger("scaff")


# ────────── helper: build bin<->scaffold DF (no temp file) ──────────
def bin_dataframe(bin_dir: Path, log: logging.Logger) -> pd.DataFrame:
    rows = []
    for fa in bin_dir.glob("*.fa"):
        bin_name = fa.stem
        for record in SeqIO.parse(fa, "fasta"):
            rows.append((record.id, bin_name))
    df = pd.DataFrame(rows, columns=["scaffold", "bin"])
    log.info(f"Collected {len(df):,} scaffold↦bin mappings from {bin_dir}")
    return df


# ────────── main pipeline ──────────
def process(scaff_tsv: Path, stb_tsv: Path, meta_csv: Path,
            bin_dir: Path, out_tsv: Path,
            len_min=1000, cov_min=5.0, br_min=0.4, log=None):

    if log is None:                         # Fallback logger
        log = logging.getLogger("scaff")

    # 1. bin map -------------------------------------------------------------
    bin_df = bin_dataframe(bin_dir, log)

    # 2. STB ----------------------------------------------------------------
    stb_df = pd.read_csv(stb_tsv, sep="\t",
                         names=["scaffold", "Sample"],  # Sample == External.ID
                         dtype={"scaffold": str, "Sample": str})

    stb_df = stb_df.merge(bin_df, on="scaffold", how="left")

    # 3. scaffold_info ------------------------------------------------------
    usecols = ["scaffold", "length", "coverage", "breadth"]
    dtypes  = {"scaffold": str,
               "length":   np.int32,
               "coverage": np.float32,
               "breadth":  np.float32}
    scaff = pd.read_csv(scaff_tsv, sep="\t", usecols=usecols, dtype=dtypes)

    log.info(f"Scaffolds loaded: {len(scaff):,}")
    scaff = scaff.query("length >= @len_min and coverage >= @cov_min and breadth >= @br_min")
    log.info(f"After quality filters: {len(scaff):,}")

    # 4. merge scaff <- stb+bin --------------------------------------------
    scaff = scaff.merge(stb_df, on="scaffold", how="inner")

    # 5. metadata  ----------------------------------------------------------
    meta_keep = [
        'External.ID', 'week_num', 'Participant ID', 'sex', 'diagnosis',
        'Height', 'Weight', 'BMI', 'fecalcal_ng_ml', 'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
        'Antibiotics', 'Immunosuppressants (e.g. oral corticosteroids)'
    ]
    meta = pd.read_csv(meta_csv, usecols=meta_keep, dtype=str)
    scaff = scaff.merge(meta, left_on="Sample", right_on="External.ID", how="left")

    # 6. summary ------------------------------------------------------------
    summary = (scaff
               .groupby("bin")[["length", "coverage", "breadth"]]
               .agg(["mean", "std", "min", "max"])
               .round(2))
    log.info("Summary per bin:\n%s", summary)

    # 7. write --------------------------------------------------------------
    scaff.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"Written {len(scaff):,} rows → {out_tsv}")


# ────────── CLI ──────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent("""\
            Merge inStrain scaffold_info with STB, bin FASTAs and sample metadata.
            Provides the same output as the original script but with less I/O
            and no intermediate files.
        """))
    ap.add_argument("--scaffold_file", required=True)
    ap.add_argument("--stb_file",      required=True)
    ap.add_argument("--metadata_file", required=True)
    ap.add_argument("--bin_dir",       required=True)
    ap.add_argument("--output_file",   required=True)
    ap.add_argument("--min_length",   type=int,   default=1000)
    ap.add_argument("--min_coverage", type=float, default=5.0)
    ap.add_argument("--min_breadth",  type=float, default=0.4)
    ap.add_argument("--log_file",     required=True)
    args = ap.parse_args()

    logger = setup_logger(args.log_file)

    process(Path(args.scaffold_file),
            Path(args.stb_file),
            Path(args.metadata_file),
            Path(args.bin_dir),
            Path(args.output_file),
            len_min=args.min_length,
            cov_min=args.min_coverage,
            br_min=args.min_breadth,
            log=logger)
