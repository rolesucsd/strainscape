#!/usr/bin/env python3
"""
Fast mapping of SNVs to Bakta genes (+nearest upstream / downstream).

This module maps mutations to gene annotations using a merge-asof
approach. For genic SNVs, it assigns the gene's
annotation. For intergenic SNVs, it merges upstream and downstream
annotations into slash-delimited fields.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import argparse
from strainscape.utils import get_logger

# Initialize logger
logger = get_logger(__name__)


def map_genes_fast(trend_path: Path, gene_path: Path, out_path: Path) -> None:
    """Map SNVs to gene/upstream/downstream annotations.

    Args:
        trend_path: TSV with columns ['scaffold','position','ref_base','new_base']
        gene_path: TSV (Bakta) with header line prefixed '#' and columns including
                   ['Sequence Id','Type','Start','Stop','Strand',
                    'Locus Tag','Gene','Product','DbXrefs']
        out_path:   Path to write annotated TSV
    """
    # Load trends
    logger.info("Loading SNV trends from %s", trend_path)
    trends = pd.read_csv(trend_path, sep="\t")
    trends = trends.rename(columns={
        'scaffold': 'Chromosome',
        'position': 'Position'
    })

    # Load gene annotations
    logger.info("Loading gene annotations from %s", gene_path)
    with open(gene_path) as gf:
        lines = gf.readlines()
    
    # Handle both formats: with # header or without
    header_lines = [ln for ln in lines if ln.startswith('#')]
    if header_lines:
        # Standard Bakta format with # header
        header_line = header_lines[-1]
        cols = header_line.lstrip('#').strip().split('\t')
        data = [ln for ln in lines if not ln.startswith('#')]
    else:
        # Format without # header - use first line as header
        header_line = lines[0]
        cols = header_line.strip().split('\t')
        data = lines[1:]
    
    # normalize column names
    cols = [c.replace(' ', '_') for c in cols]
    from io import StringIO
    genes = pd.read_csv(StringIO(''.join(data)), sep='\t', names=cols)
    
    # Use 'scaffold' column if available (full scaffold name), otherwise use 'Sequence_Id'
    if 'scaffold' in genes.columns:
        genes = genes.rename(columns={'scaffold': 'Chromosome'})
        logger.info("Using 'scaffold' column for chromosome matching")
    elif 'Sequence_Id' in genes.columns:
        genes = genes.rename(columns={'Sequence_Id': 'Chromosome'})
        logger.info("Using 'Sequence_Id' column for chromosome matching")
    else:
        raise ValueError("Gene file must have either 'scaffold' or 'Sequence Id' column")

    # ensure numeric ordering
    trends['Position'] = trends['Position'].astype(int)
    genes['Start'] = genes['Start'].astype(int)
    genes['Stop'] = genes['Stop'].astype(int)

    # annotation fields to merge
    annot_cols = ['Type','Start','Stop','Strand','Locus_Tag','Gene','Product','DbXrefs']

    all_results = []
    # process per chromosome
    for chrom, snvs in trends.groupby('Chromosome', sort=False):
        chr_genes = genes.loc[genes['Chromosome'] == chrom]
        if chr_genes.empty:
            logger.warning("No genes for %s; skipping", chrom)
            continue
        # sort genes by Start
        chr_genes = chr_genes.sort_values('Start').reset_index(drop=True)
        # sort SNVs by Position
        snvs = snvs.sort_values('Position').reset_index(drop=True)

        # prepare upstream / downstream tables
        g_up = chr_genes.rename(columns={c: c + '_up' for c in annot_cols})
        g_dn = chr_genes.rename(columns={c: c + '_dn' for c in annot_cols})

        # nearest upstream: Start_up <= Position
        up = pd.merge_asof(
            snvs, g_up,
            left_on='Position', right_on='Start_up',
            by='Chromosome',
            direction='backward'
        )
        # nearest downstream: Start_dn > Position
        dn = pd.merge_asof(
            snvs, g_dn,
            left_on='Position', right_on='Start_dn',
            by='Chromosome',
            direction='forward'
        )
        # combine
        merged = up.join(dn[[col for col in dn.columns if col.endswith('_dn')]])

        # detect genic: Position between Start_up and Stop_up
        in_gene = merged['Position'].between(
            merged['Start_up'], merged['Stop_up']
        )
        merged['gene_type'] = np.where(in_gene, 'genic', 'intergenic')

        # for each annotation, choose or merge
        for col in annot_cols:
            up_col = col + '_up'
            dn_col = col + '_dn'
            merged[col] = np.where(
                in_gene,
                merged[up_col].astype(str),
                merged[up_col].astype(str) + '/' + merged[dn_col].astype(str)
            )

        # Keep all original columns plus the new annotation columns
        # Preserve depth columns (depth_mean, depth_std, depth_min, depth_max, position_coverage, etc.)
        cols_out = [col for col in snvs.columns] + ['gene_type'] + annot_cols
        # Only include columns that actually exist
        cols_out = [col for col in cols_out if col in merged.columns]
        all_results.append(merged[cols_out])

    # concatenate and write
    if all_results:
        df_out = pd.concat(all_results, ignore_index=True)
    else:
        # empty structure with all original columns plus annotations
        cols_out = [col for col in trends.columns] + ['gene_type'] + annot_cols
        df_out = pd.DataFrame(columns=cols_out)
        logger.warning("No SNV annotations; writing empty file %s", out_path)

    df_out.to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote %d annotations to %s", len(df_out), out_path)


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description='Annotate SNVs as genic vs intergenic.'
    )
    p.add_argument('--trend_file', required=True,
                   help='TSV of SNVs [scaffold pos ref_base new_base]')
    p.add_argument('--gene_file', required=True,
                   help='Bakta TSV with header # line and gene columns')
    p.add_argument('--output_file', required=True,
                   help='Path to output annotated TSV')
    args = p.parse_args()
    map_genes_fast(
        Path(args.trend_file),
        Path(args.gene_file),
        Path(args.output_file)
    )