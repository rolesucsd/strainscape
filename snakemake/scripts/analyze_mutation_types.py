#!/usr/bin/env python3
"""
Annotate SNVs as Silent / Missense / Nonsense
(keeps Renee's original logic & quirky +1 index, just faster).
"""

from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import functools, logging, argparse, sys

# ───────── logging ─────────
def setup_logger(logf=None):
    fmt = "%(asctime)s %(levelname)s  %(message)s"
    hdlr = logging.FileHandler(logf) if logf else logging.StreamHandler(sys.stderr)
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[hdlr])
    return logging.getLogger("mut_annot")

log = setup_logger()

# ───────── core annotator (same biology, cached) ─────────
def annotate_row(row, seqs, gene_cache):
    """
    Return (mutation_type, coding_status) for one row.

    Keeps Renee's original +1 position offset so that base calls align
    with her tables.  All other logic is unchanged.
    """
    # intergenic fast-exit
    if row.get("gene_type") == "intergenic":
        return "Intergenic", "Non-Coding"

    start  = int(row["Start"])
    stop   = int(row["Stop"])
    strand = row["Strand"]
    pos1   = int(row["Position"])          # 1-based
    ref_b  = row["ref_base"].upper()
    alt_b  = row["new_base"].upper()
    chrom  = row["Chromosome"]

    seq = seqs.get(chrom)
    if seq is None:
        return "Seq_not_found", "Error"

    key = (chrom, start, stop)
    gene_fwd = gene_cache.get(key)
    if gene_fwd is None:
        gene_fwd = seq[start-1:stop]       # forward string
        gene_cache[key] = gene_fwd

    # **original offset (+1) kept here**
    idx = pos1 - start + 1                 # 0-based in gene_fwd
    if idx < 0 or idx >= len(gene_fwd):
        return "Pos_out_of_bounds", "Error"

    if gene_fwd[idx].upper() != ref_b:
        return f"Reference base mismatch (expected {ref_b}, got {gene_fwd[idx]})", "Error"

    # build mutated codon (still forward strand)
    mut_fwd = gene_fwd[:idx] + alt_b + gene_fwd[idx+1:]
    codon0  = (idx // 3) * 3
    ref_cdn = gene_fwd[codon0:codon0+3]
    alt_cdn = mut_fwd[codon0:codon0+3]

    if strand == '-':
        ref_cdn = str(Seq(ref_cdn).reverse_complement())
        alt_cdn = str(Seq(alt_cdn).reverse_complement())

    ref_aa = str(Seq(ref_cdn).translate())
    alt_aa = str(Seq(alt_cdn).translate())

    if ref_aa == alt_aa:
        return "Silent", "Coding"
    if alt_aa == '*' and ref_aa != '*':
        return "Nonsense", "Coding"
    return "Missense", "Coding"


# ───────── main wrapper ─────────
def run(muts_tsv: Path, ref_fa: Path, out_tsv: Path):

    # read reference
    log.info("Loading FASTA …")
    seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(ref_fa, "fasta")}
    if not seqs:
        raise RuntimeError("reference FASTA empty")

    # read table
    log.info("Loading mutations …")
    muts = pd.read_csv(muts_tsv, sep="\t")

    # annotate
    cache = {}
    fn = functools.partial(annotate_row, seqs=seqs, gene_cache=cache)
    log.info("Annotating …")
    mut_type, coding = zip(*muts.apply(fn, axis=1))
    muts["Mutation_Type"] = mut_type
    muts["Coding_Status"] = coding

    muts.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"Saved → {out_tsv}")


# ───────── CLI ─────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Annotate SNVs (silent / missense / nonsense)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mutation_file",  required=True)
    ap.add_argument("--reference_file", required=True)
    ap.add_argument("--output_file",    required=True)
    ap.add_argument("--log_file")
    args = ap.parse_args()

    if args.log_file:
        log = setup_logger(args.log_file)

    run(Path(args.mutation_file),
        Path(args.reference_file),
        Path(args.output_file))
