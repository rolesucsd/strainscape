#!/usr/bin/env python3
"""
Build a Prodigal gene-coordinate map from contigs.proteins.faa.

Input headers look like:
>bc2203_MaxBin_bin.14.fa|s52.ctg000055c_1 # 107 # 1387 # 1 # ID=...

Output columns:
- scaffold_full  (e.g., bc2203_MaxBin_bin.14.fa|s52.ctg000055c)
- bin            (left of '|')
- scaffold       (right of '|')
- gene_idx       (int, e.g., 1)
- start, stop    (1-based inclusive)
- strand         (as given by Prodigal: 1 or -1)
- gene_id_full   (= scaffold_full + "_" + gene_idx)

File type is inferred from --out: .parquet / .feather / .tsv
"""

from pathlib import Path
import argparse, re
import pandas as pd

def write_table(df: pd.DataFrame, out: Path):
    p = str(out)
    if p.endswith(".parquet"):
        df.to_parquet(p, index=False)
    elif p.endswith(".feather"):
        df.to_feather(p)
    else:
        df.to_csv(p, sep="\t", index=False)

def main(faa: Path, out: Path):
    rows = []
    with open(faa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # keep only the first token (before any spaces)
            first = line[1:].strip().split()[0]
            # split the '#'-delimited metadata
            parts = [p.strip() for p in line[1:].split("#")]
            # contig_geneIdx is in the first part
            m = re.match(r"(.+)_([0-9]+)$", first)
            if not m:
                continue
            scaffold_full = m.group(1)
            gene_idx      = int(m.group(2))
            # Prodigal puts start, stop, strand in the next 3 fields
            start = int(parts[1].strip())
            stop  = int(parts[2].strip())
            strand = parts[3].strip()  # '1' or '-1'
            if start > stop:  # just in case
                start, stop = stop, start
            # bin|scaffold split
            if "|" in scaffold_full:
                bin_, scaf = scaffold_full.split("|", 1)
            else:
                bin_, scaf = scaffold_full, scaffold_full
            
            gene_id_full = f"{scaffold_full}_{gene_idx}"
            
            rows.append({
                'scaffold_full': scaffold_full,
                'bin': bin_,
                'scaffold': scaf,
                'gene_idx': gene_idx,
                'start': start,
                'stop': stop,
                'strand': strand,
                'gene_id_full': gene_id_full
            })
    
    df = pd.DataFrame(rows)
    write_table(df, out)
    print(f"Parsed {len(df)} genes from {faa}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Prodigal protein FASTA headers")
    parser.add_argument("--faa", type=Path, required=True, help="Input protein FASTA file")
    parser.add_argument("--out", type=Path, required=True, help="Output file (.tsv, .parquet, or .feather)")
    
    args = parser.parse_args()
    main(args.faa, args.out)
