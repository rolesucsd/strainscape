#!/usr/bin/env python3
"""
Annotate SNVs as Silent / Missense / Nonsense
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

def load_sequences(fasta_file: Path) -> dict:
    """Load sequences from a FASTA file into a dictionary."""
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_file, "fasta")}

def get_mutation_type(row: pd.Series, sequences: dict) -> dict:
    """
    Determine the type of mutation for a single row.
    
    Args:
        row: DataFrame row containing mutation information
        sequences: Dictionary of chromosome sequences
        
    Returns:
        Dictionary with mutation type information
    """
    # intergenic fast-exit
    if row.get("coding") == "intergenic":
        return {"Mutation_Type": "Silent", "Coding_Status": "Non-Coding"}

    start = int(row["Matched_Start"])
    stop = int(row["Matched_Stop"])
    strand = row["Matched_Strand"]
    pos1 = int(row["Position"])  # 1-based
    ref_b = row["ref_base"].upper()
    alt_b = row["new_base"].upper()
    chrom = row["Chromosome"]

    seq = sequences.get(chrom)
    if seq is None:
        return {"Mutation_Type": "Silent", "Coding_Status": "Error"}

    # Get gene sequence
    gene_fwd = seq[start-1:stop]  # forward string

    # **original offset (+1) kept here**
    idx = pos1 - start + 1  # 0-based in gene_fwd
    if idx < 0 or idx >= len(gene_fwd):
        return {"Mutation_Type": "Silent", "Coding_Status": "Error"}

    if gene_fwd[idx].upper() != ref_b:
        return {"Mutation_Type": "Silent", "Coding_Status": "Error"}

    # build mutated codon (still forward strand)
    mut_fwd = gene_fwd[:idx] + alt_b + gene_fwd[idx+1:]
    codon0 = (idx // 3) * 3
    ref_cdn = gene_fwd[codon0:codon0+3]
    alt_cdn = mut_fwd[codon0:codon0+3]

    if strand == '-':
        ref_cdn = str(Seq(ref_cdn).reverse_complement())
        alt_cdn = str(Seq(alt_cdn).reverse_complement())

    ref_aa = str(Seq(ref_cdn).translate())
    alt_aa = str(Seq(alt_cdn).translate())

    if ref_aa == alt_aa:
        return {"Mutation_Type": "Silent", "Coding_Status": "Coding"}
    if alt_aa == '*' and ref_aa != '*':
        return {"Mutation_Type": "Nonsense", "Coding_Status": "Coding"}
    return {"Mutation_Type": "Missense", "Coding_Status": "Coding"}

def analyze_mutation_types(mutations: pd.DataFrame, sequences: dict) -> pd.DataFrame:
    """Analyze mutation types for a DataFrame of mutations.
    
    Args:
        mutations: DataFrame containing mutation data
        sequences: Dictionary mapping chromosome names to sequences
        
    Returns:
        DataFrame with mutation type annotations
    """
    # Create a copy to avoid modifying the input
    result = mutations.copy()
    if result.empty:
        result['Mutation_Type'] = []
        result['Coding_Status'] = []
        return result
    # Apply mutation type analysis
    fn = functools.partial(get_mutation_type, sequences=sequences)
    mut_type, coding = zip(*result.apply(fn, axis=1))
    result["Mutation_Type"] = mut_type
    result["Coding_Status"] = coding
    return result

# ───────── main wrapper ─────────
def run(muts_tsv: Path, ref_fa: Path, out_tsv: Path):

    # read reference
    log.info("Loading FASTA …")
    seqs = load_sequences(ref_fa)
    if not seqs:
        raise RuntimeError("reference FASTA empty")

    # read table
    log.info("Loading mutations …")
    muts = pd.read_csv(muts_tsv, sep="\t")

    # annotate
    cache = {}
    fn = functools.partial(get_mutation_type, sequences=seqs)
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
