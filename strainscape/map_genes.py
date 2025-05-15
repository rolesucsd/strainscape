#!/usr/bin/env python3
"""
Fast mapping of SNVs to Bakta genes (+nearest upstream / downstream).

This module provides functions to map mutations to genes and their surrounding
regions using Bakta gene annotations. It identifies whether mutations are
genic or intergenic and finds the nearest upstream and downstream genes.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import logging, argparse, sys
from typing import Dict, List, Optional, Tuple
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

# ──────────────────────── core function ─────────────────────
def map_genes_fast(mut_file: Path,
                   trend_file: Path,
                   gene_file: Path,
                   out_file: Path) -> None:
    """Map mutations to genes and their surrounding regions.
    
    Args:
        mut_file: Path to TSV file containing mutation data with columns:
            - Chromosome: str, chromosome/scaffold name
            - Position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
            - A/C/G/T: int, base counts
            - position_coverage: int, total coverage
        trend_file: Path to TSV file containing trend data with columns:
            - scaffold: str, chromosome/scaffold name
            - position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
        gene_file: Path to Bakta gene annotation file with columns:
            - Sequence Id: str, chromosome/scaffold name
            - Type: str, feature type (e.g. "CDS")
            - Start: int, gene start position
            - Stop: int, gene stop position
            - Strand: str, "+" or "-"
            - Locus Tag: str, gene identifier
            - Gene: str, gene name
            - Product: str, gene product
        out_file: Path to output TSV file.
        
    Returns:
        None. Writes results to out_file with columns:
            - Chromosome: str, chromosome/scaffold name
            - Position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
            - gene_type: str, "genic" or "intergenic"
            - Type: str, gene type
            - Start: int, gene start position
            - Stop: int, gene stop position
            - Strand: str, "+" or "-"
            - Locus Tag: str, gene identifier
            - Gene: str, gene name
            - Product: str, gene product
            - Start_up: int, upstream gene start
            - Stop_up: int, upstream gene stop
            - Strand_up: str, upstream gene strand
            - Locus Tag_up: str, upstream gene identifier
            - Gene_up: str, upstream gene name
            - Product_up: str, upstream gene product
            - Start_dn: int, downstream gene start
            - Stop_dn: int, downstream gene stop
            - Strand_dn: str, downstream gene strand
            - Locus Tag_dn: str, downstream gene identifier
            - Gene_dn: str, downstream gene name
            - Product_dn: str, downstream gene product
    """
    logger.info("Loading data …")
    muts   = pd.read_csv(mut_file,   sep="\t")
    trends = pd.read_csv(trend_file, sep="\t")
    
    # Read gene data with proper header handling
    with open(gene_file, 'r') as f:
        lines = f.readlines()
    header_line = [line for line in lines if line.startswith('#')][-1]
    header = header_line.lstrip('#').strip().split('\t')
    data_lines = [line for line in lines if not line.startswith('#')]
    from io import StringIO
    data_str = ''.join(data_lines)
    genes = pd.read_csv(StringIO(data_str), sep='\t', names=header)

    # harmonise column names once
    trends = trends.rename(columns={"scaffold": "Chromosome",
                                    "position": "Position"})
    muts   = pd.merge(muts, trends,
                      on=["Chromosome", "Position", "ref_base", "new_base"],
                      how="left")

    # ensure numeric
    muts["Position"]   = muts["Position"].astype(int)
    genes["Start"]     = genes["Start"].astype(int)
    genes["Stop"]      = genes["Stop"].astype(int)

    # containers
    out_frames = []

    # ───────── process chromosome by chromosome ─────────
    for chrom, mut_chunk in muts.groupby("Chromosome", sort=False):
        g = genes.loc[genes["Sequence Id"] == chrom]
        if g.empty:
            logger.warning(f"No genes for {chrom}; skipping")
            continue

        # sort once for merge_asof
        g = g.sort_values("Start").reset_index(drop=True)
        m = mut_chunk.sort_values("Position").reset_index(drop=True)

        # nearest upstream (gene start <= pos)
        up  = pd.merge_asof(m, g,
                            left_on="Position", right_on="Start",
                            direction="backward",
                            suffixes=("", "_up"))

        # nearest downstream (gene start > pos)
        dn  = pd.merge_asof(m, g,
                            left_on="Position", right_on="Start",
                            direction="forward",
                            suffixes=("", "_dn"))

        # combine
        merged = up.join(dn.filter(regex="_dn$"))

        # identify whether the mutation is inside the upstream gene
        in_gene = merged["Position"].between(merged["Start"], merged["Stop"])
        merged.loc[in_gene, "gene_type"] = merged.loc[in_gene, "Type"]   # genic
        merged.loc[~in_gene, "gene_type"] = "intergenic"

        out_frames.append(merged)

    # ───────── final tidy-up & save ─────────
    result = (pd.concat(out_frames, ignore_index=True)
                .rename(columns={"Type": "mutation_type"}))

    result.to_csv(out_file, sep="\t", index=False)
    logger.info(f"Written {len(result):,} rows → {out_file}")


# ────────────────────────── CLI ─────────────────────────────
if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--mutation_file", required=True)
    p.add_argument("--trend_file",    required=True)
    p.add_argument("--gene_file",     required=True)
    p.add_argument("--output_file",   required=True)
    args = p.parse_args()

    map_genes_fast(Path(args.mutation_file),
                   Path(args.trend_file),
                   Path(args.gene_file),
                   Path(args.output_file))