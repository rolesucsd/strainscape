#!/usr/bin/env python3
"""
SNV → gene mapping (refactor):
- Reads a prebuilt Prodigal map (no FASTA parsing here)
- Expands only intergenic SNVs to two rows (up & dn)
- Merges clusters (gene_id_full → cluster_rep) and EggNOG (cluster_rep → anno)
- Adds Bakta nearest upstream/downstream annotations, selecting the flank-specific set per row
- Prints debug row counts at each join

Inputs
------
--trend_file       TSV/Parquet/Feather with columns:
                   ['scaffold','position','ref_base','new_base'] where scaffold="bin|contig"
--prodigal_map     TSV/Parquet/Feather with columns:
                   ['scaffold_full','start','stop','gene_id_full'] (others OK but not required)
--clusters_file    2-col TSV/Parquet/Feather:
                   [cluster_rep, gene_id_full]  (deduped to 1:1 inside script)
--eggnog_file      TSV/Parquet/Feather with column:
                   ['query'] (unique per representative; deduped to 1:1 inside script)
--gene_file        Bakta TSV (may start with '#'), includes Sequence Id / Start / Stop / ...
--output_file      Path; type inferred from extension (.parquet/.feather/.tsv)

Outputs
-------
One row per genic SNV; two rows per intergenic SNV (up & dn).
Keeps: bin, scaffold (post-split), scaffold_full, gene/flanks, cluster_rep, EggNOG,
       Bakta_* (flank-selected).

Debug prints show row counts before/after each join.
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
from io import StringIO

# ---------------- I/O helpers (minimal, fast) ----------------

def read_table(p: Path, sep="\t"):
    p = str(p)
    if p.endswith(".parquet"):
        return pd.read_parquet(p)
    if p.endswith(".feather"):
        return pd.read_feather(p)
    try:
        return pd.read_csv(p, sep=sep, engine="pyarrow")
    except Exception:
        return pd.read_csv(p, sep=sep)

def write_table(df: pd.DataFrame, p: Path):
    p = str(p)
    if p.endswith(".parquet"):
        df.to_parquet(p, index=False)
    elif p.endswith(".feather"):
        df.to_feather(p)
    else:
        df.to_csv(p, sep="\t", index=False)

def read_bakta_genes(path: Path) -> pd.DataFrame:
    # supports Bakta header lines starting with '#'
    with open(path) as f:
        lines = f.readlines()
    if lines and lines[0].startswith('#'):
        header = [ln for ln in lines if ln.startswith('#')][-1].lstrip('#').strip()
        cols = [c.replace(' ', '_') for c in header.split('\t')]
        data = ''.join(ln for ln in lines if not ln.startswith('#'))
        genes = pd.read_csv(StringIO(data), sep='\t', names=cols)
    else:
        genes = read_table(path)
        genes.columns = [c.replace(' ', '_') for c in genes.columns]
    return genes

# ---------------- Core ----------------

def annotate(trend_file: Path,
             prodigal_map: Path,
             clusters_file: Path,
             eggnog_file: Path,
             gene_file: Path,
             output_file: Path):

    # Trends
    trends = read_table(trend_file)
    print(f"[DEBUG] trends read: rows={len(trends)}")
    trends = trends.rename(columns={'scaffold': 'scaffold_full'})
    # Remove .fa extension from scaffold names to match Bakta format
    trends['scaffold_full'] = trends['scaffold_full'].str.replace('.fa|', '|', regex=False)
    trends[['bin','scaffold']] = trends['scaffold_full'].str.split('|', n=1, expand=True)

    # Prodigal gene map (prebuilt)
    gmap = read_table(prodigal_map)
    print(f"[DEBUG] prodigal_map read: rows={len(gmap)}")
    # Remove .fa extension from scaffold names to match trends format
    gmap['scaffold_full'] = gmap['scaffold_full'].str.replace('.fa|', '|', regex=False)
    # Ensure required cols exist
    gmap = gmap[['scaffold_full','start','stop','gene_id_full']].copy()
    # group per scaffold for fast lookups
    g_by = {s: g.sort_values('start').reset_index(drop=True)
            for s, g in gmap.groupby('scaffold_full', sort=False)}

    # Genic/flank assignment using Prodigal map
    out_parts = []
    for scaf, sn in trends.groupby('scaffold_full', sort=False):
        sn = sn.copy()
        g = g_by.get(scaf)
        if g is None or g.empty:
            sn['gene_id_full'] = pd.NA
            sn['gene_id_full_up'] = pd.NA
            sn['gene_id_full_dn'] = pd.NA
            sn['intergenic'] = 'upstream'  # no genes present; arbitrary
            out_parts.append(sn)
            continue

        starts = g['start'].to_numpy()
        stops  = g['stop'].to_numpy()
        gids   = g['gene_id_full'].to_numpy()

        pos = sn['position'].to_numpy()
        j = np.searchsorted(starts, pos, side='right')  # first gene start > pos (downstream idx)
        i = j - 1                                       # upstream idx (<= pos)

        # genic if starts[i] <= pos <= stops[i]
        in_range = (
            (i >= 0) & (i < len(starts)) &
            (starts[np.clip(i,0,len(starts)-1)] <= pos) &
            (pos <= stops[np.clip(i,0,len(stops)-1)])
        )
        gene_id_full = np.where(
            in_range, gids[np.clip(i,0,len(gids)-1)], None
        )

        gene_up = np.where((i >= 0) & (i < len(gids)), gids[np.clip(i,0,len(gids)-1)], None)
        gene_dn = np.where((j >= 0) & (j < len(gids)), gids[np.clip(j,0,len(gids)-1)], None)

        # nearest side for intergenic flag
        d_up = np.where((i >= 0) & (i < len(stops)), np.abs(pos - stops[np.clip(i,0,len(stops)-1)]), np.inf)
        d_dn = np.where((j >= 0) & (j < len(starts)), np.abs(starts[np.clip(j,0,len(starts)-1)] - pos), np.inf)
        intergenic = np.where(in_range, 'genic', np.where(d_dn <= d_up, 'upstream', 'downstream'))

        sn['gene_id_full'] = gene_id_full
        sn['gene_id_full_up'] = gene_up
        sn['gene_id_full_dn'] = gene_dn
        sn['intergenic'] = intergenic
        out_parts.append(sn)

    snv0 = pd.concat(out_parts, ignore_index=True)
    print(f"[DEBUG] after Prodigal mapping: rows={len(snv0)} "
          f"(genic={snv0['intergenic'].eq('genic').sum()}, intergenic={(snv0['intergenic']!='genic').sum()})")

    # Expand ONLY intergenic rows to two rows (up & dn). Genic stays single.
    genic = snv0[snv0['intergenic'].eq('genic')].copy()
    genic['target_gene'] = genic['gene_id_full']
    genic['flank'] = 'self'

    inter = snv0[snv0['intergenic'].ne('genic')].copy()

    up = inter.copy()
    up['target_gene'] = up['gene_id_full_up']
    up['flank'] = 'up'
    up['intergenic'] = 'upstream'

    dn = inter.copy()
    dn['target_gene'] = dn['gene_id_full_dn']
    dn['flank'] = 'dn'
    dn['intergenic'] = 'downstream'

    snv = pd.concat([genic, up, dn], ignore_index=True)
    snv = snv[snv['target_gene'].notna()].reset_index(drop=True)

    print(f"[DEBUG] after intergenic expansion: rows={len(snv)} "
          f"(expected ≈ genic + 2*intergenic_with_flanks)")

    # ------ Merge CLUSTERS (many_to_one) ------
    clusters = read_table(clusters_file)
    if clusters.shape[1] == 2:
        clusters.columns = ['cluster_rep','gene_id_full']  # <- matches your working input
    dup_before = clusters.duplicated(subset=['gene_id_full']).sum()
    clusters = clusters.drop_duplicates(subset=['gene_id_full'], keep='first')
    dup_after = clusters.duplicated(subset=['gene_id_full']).sum()
    print(f"[DEBUG] clusters read: rows={len(clusters)} (dups_by_gene_id_full before={dup_before}, after={dup_after})")

    rows_before = len(snv)
    snv = snv.merge(clusters.rename(columns={'gene_id_full':'target_gene'}),
                    on='target_gene', how='left', validate='many_to_one')
    print(f"[DEBUG] after merge clusters: rows={len(snv)} (delta={len(snv)-rows_before})")

    # ------ Merge EGGNOG (many_to_one) ------
    egg = read_table(eggnog_file)
    if 'query' not in egg.columns:
        raise ValueError("EggNOG file must have a 'query' column.")
    dup_egg_before = egg.duplicated(subset=['query']).sum()
    egg = egg.drop_duplicates(subset=['query'], keep='first')
    dup_egg_after = egg.duplicated(subset=['query']).sum()
    reps = snv['cluster_rep'].dropna().astype(str).unique()
    egg = egg[egg['query'].astype(str).isin(reps)]
    print(f"[DEBUG] eggnog read: rows={len(egg)} "
          f"(dups_by_query before={dup_egg_before}, after={dup_egg_after}, "
          f"prefiltered_to_present_reps={len(egg)})")

    rows_before = len(snv)
    snv = snv.merge(egg, left_on='cluster_rep', right_on='query',
                    how='left', validate='many_to_one')
    print(f"[DEBUG] after merge eggnog: rows={len(snv)} (delta={len(snv)-rows_before})")

    # ------ Bakta nearest upstream/downstream (by post-split scaffold) ------
    bak = read_bakta_genes(gene_file)
    need = ['scaffold_full','Type','Start','Stop','Strand','Locus_Tag','Gene','Product','DbXrefs']
    for c in need:
        if c not in bak.columns:
            bak[c] = pd.NA
    bak = bak[need].copy()
    bak['Start'] = pd.to_numeric(bak['Start'], errors='coerce').astype(np.int64)
    bak['Stop']  = pd.to_numeric(bak['Stop'],  errors='coerce').astype(np.int64)

    print(f"[DEBUG] bakta read: rows={len(bak)} (unique scaffolds={bak['scaffold_full'].nunique()})")

    parts = []
    for scaf, part in snv.groupby('scaffold_full', sort=False):
        g = bak[bak['scaffold_full'] == scaf].sort_values('Start').reset_index(drop=True)
        if g.empty:
            part = part.copy()
            for c in need[1:]:
                part[c + '_up'] = pd.NA
                part[c + '_dn'] = pd.NA
            parts.append(part)
            continue

        g_up = g.rename(columns={c: f"{c}_up" for c in need[1:]})
        g_dn = g.rename(columns={c: f"{c}_dn" for c in need[1:]})
        part = part.sort_values('position')
        up = pd.merge_asof(part, g_up, left_on='position', right_on='Start_up',
                           by='scaffold_full', direction='backward')
        dn = pd.merge_asof(part, g_dn, left_on='position', right_on='Start_dn',
                           by='scaffold_full', direction='forward')
        # join only the downstream-specific columns to avoid name clashes
        dn_cols = [c for c in dn.columns if c.endswith('_dn') and c not in up.columns]
        merged = up.join(dn[dn_cols]) if dn_cols else up
        parts.append(merged)

    rows_before = len(snv)
    snv = pd.concat(parts, ignore_index=True)
    print(f"[DEBUG] after Bakta nearest up/dn: rows={len(snv)} (delta={len(snv)-rows_before})")

    # ---- Select the flank-appropriate Bakta fields into a single set of columns ----
    def pick(col):
        return np.where(snv['flank'].eq('up'), snv[f'{col}_up'],
                        np.where(snv['flank'].eq('dn'), snv[f'{col}_dn'], snv[f'{col}_up']))

    for col in ['Type','Strand','Locus_Tag','Gene','Product','DbXrefs','Start','Stop']:
        snv[f'Bakta_{col}'] = pick(col)

    # Drop the intermediate *_up / *_dn Bakta columns to keep only flank-selected Bakta_*
    bakta_tmp_cols = [c for c in snv.columns if c.endswith('_up') or c.endswith('_dn')]
    snv = snv.drop(columns=bakta_tmp_cols, errors='ignore')
    print(f"[DEBUG] dropped Bakta *_up/*_dn columns: n={len(bakta_tmp_cols)}")

    # Final column order
    base = ['scaffold_full','bin','scaffold','position','ref_base','new_base',
            'gene_id_full','gene_id_full_up','gene_id_full_dn','intergenic',
            'target_gene','cluster_rep']
    base = [c for c in base if c in snv.columns]
    # Keep EggNOG columns (already merged once)
    other = [c for c in snv.columns if c not in set(base + ['query'] + [c for c in snv.columns if c.startswith('Bakta_')])]
    bakta_sel = [c for c in snv.columns if c.startswith('Bakta_')]
    final_cols = base + other + ['query'] + bakta_sel

    write_table(snv[final_cols], output_file)
    print(f"[DEBUG] wrote: rows={len(snv)} to {output_file}")

# ---------------- CLI ----------------

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="SNV→gene mapping using prebuilt Prodigal map, clusters, EggNOG, and Bakta (flank-selected).")
    ap.add_argument("--trend_file",    required=True)
    ap.add_argument("--prodigal_map",  required=True, help="Output of build_prodigal_map.py")
    ap.add_argument("--clusters_file", required=True)
    ap.add_argument("--eggnog_file",   required=True)
    ap.add_argument("--gene_file",     required=True, help="Bakta TSV (may have # header lines)")
    ap.add_argument("--output_file",   required=True)
    args = ap.parse_args()

    annotate(
        Path(args.trend_file),
        Path(args.prodigal_map),
        Path(args.clusters_file),
        Path(args.eggnog_file),
        Path(args.gene_file),
        Path(args.output_file),
    )
