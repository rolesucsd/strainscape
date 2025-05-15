#!/usr/bin/env python3
"""
Analyze mutation types and patterns.

This module provides functions to analyze mutation types and patterns in
genomic data. It processes mutation data to identify different types of
mutations (e.g., transitions, transversions) and their frequencies.

Inputs:
  - mutations.tsv: Filtered mutation data
Outputs:
  - mutation_types.tsv: Summary of mutation types and frequencies
"""

from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import functools, logging, argparse, sys
from typing import Dict, Tuple, Optional
from strainscape.utils import setup_logging, get_logger

# Get module logger
logger = get_logger(__name__)

def load_sequences(fasta_file: Path) -> Dict[str, str]:
    """Load sequences from a FASTA file into a dictionary.
    
    Args:
        fasta_file: Path to FASTA file containing reference sequences.
        
    Returns:
        Dictionary mapping sequence IDs to uppercase DNA sequences.
        
    Raises:
        RuntimeError: If no sequences are found in the FASTA file.
    """
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_file, "fasta")}

def get_mutation_type(row: pd.Series, sequences: Dict[str, str], translation_table: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """Determine the type of mutation for a single row.
    
    Args:
        row: DataFrame row containing mutation information with columns:
            - Chromosome: str, chromosome/scaffold name
            - Position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
            - coding: str, "genic" or "intergenic"
            - Matched_Start: int, gene start position
            - Matched_Stop: int, gene stop position
            - Matched_Strand: str, "+" or "-"
        sequences: Dictionary mapping chromosome names to DNA sequences.
        translation_table: Optional codon->amino acid translation table.
    Returns:
        Tuple (Mutation_Type, Coding_Status)
    """
    # intergenic fast-exit
    if row.get("coding") == "intergenic":
        return ("Intergenic", "Non-Coding")

    start = int(row["Matched_Start"])
    stop = int(row["Matched_Stop"])
    strand = row["Matched_Strand"]
    pos1 = int(row["Position"])  # 1-based
    ref_b = row["ref_base"].upper()
    alt_b = row["new_base"].upper()
    chrom = row["Chromosome"]

    seq = sequences.get(chrom)
    if seq is None:
        return ("Silent", "Error")

    gene_fwd = seq[start-1:stop]  # forward string
    idx = pos1 - start + 1  # 0-based in gene_fwd
    if idx < 0 or idx >= len(gene_fwd):
        return ("Silent", "Error")
    if gene_fwd[idx].upper() != ref_b:
        return ("Silent", "Error")

    mut_fwd = gene_fwd[:idx] + alt_b + gene_fwd[idx+1:]
    codon0 = (idx // 3) * 3
    ref_cdn = gene_fwd[codon0:codon0+3]
    alt_cdn = mut_fwd[codon0:codon0+3]

    if translation_table is None:
        from Bio.Data import CodonTable
        translation_table = CodonTable.unambiguous_dna_by_id[1].forward_table
        translation_table = translation_table.copy()
        translation_table["TAA"] = translation_table["TAG"] = translation_table["TGA"] = "*"

    def translate_codon(codon: str) -> str:
        codon = codon.upper().replace("U", "T")
        return translation_table.get(codon, "X")

    if strand == '-':
        ref_cdn = str(Seq(ref_cdn).reverse_complement())
        alt_cdn = str(Seq(alt_cdn).reverse_complement())

    ref_aa = translate_codon(ref_cdn)
    alt_aa = translate_codon(alt_cdn)

    if ref_aa == alt_aa:
        return ("Silent", "Coding")
    if alt_aa == '*' and ref_aa != '*':
        return ("Nonsense", "Coding")
    return ("Missense", "Coding")

def analyze_mutation_types(mutations: pd.DataFrame, sequences: Dict[str, str]) -> pd.DataFrame:
    """Analyze mutation types for a DataFrame of mutations.
    
    Args:
        mutations: DataFrame containing mutation data with columns:
            - Chromosome: str, chromosome/scaffold name
            - Position: int, 1-based position
            - ref_base: str, reference base
            - new_base: str, alternate base
            - coding: str, "genic" or "intergenic"
            - Matched_Start: int, gene start position
            - Matched_Stop: int, gene stop position
            - Matched_Strand: str, "+" or "-"
        sequences: Dictionary mapping chromosome names to DNA sequences.
        
    Returns:
        DataFrame with added columns:
            - Mutation_Type: str, one of "Silent", "Missense", "Nonsense"
            - Coding_Status: str, one of "Coding", "Non-Coding", "Error"
    """
    result = mutations.copy()
    if result.empty:
        result['Mutation_Type'] = []
        result['Coding_Status'] = []
        return result
    from Bio.Data import CodonTable
    translation_table = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
    translation_table["TAA"] = translation_table["TAG"] = translation_table["TGA"] = "*"
    fn = functools.partial(get_mutation_type, sequences=sequences, translation_table=translation_table)
    mut_type, coding = zip(*result.apply(fn, axis=1))
    result["Mutation_Type"] = mut_type
    result["Coding_Status"] = coding
    return result

def run(muts_tsv: Path, ref_fa: Path, out_tsv: Path) -> None:
    """Run mutation type analysis on a TSV file of mutations.
    
    Args:
        muts_tsv: Path to TSV file containing mutation data.
        ref_fa: Path to FASTA file containing reference sequences.
        out_tsv: Path to output TSV file.
        
    Raises:
        RuntimeError: If reference FASTA is empty.
    """
    logger.info("Loading FASTA …")
    seqs = load_sequences(ref_fa)
    if not seqs:
        raise RuntimeError("reference FASTA empty")

    logger.info("Loading mutations …")
    muts = pd.read_csv(muts_tsv, sep="\t")

    from Bio.Data import CodonTable
    translation_table = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
    translation_table["TAA"] = translation_table["TAG"] = translation_table["TGA"] = "*"
    fn = functools.partial(get_mutation_type, sequences=seqs, translation_table=translation_table)
    logger.info("Annotating …")
    mut_type, coding = zip(*muts.apply(fn, axis=1))
    muts["Mutation_Type"] = mut_type
    muts["Coding_Status"] = coding

    muts.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Saved → {out_tsv}")


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
        setup_logging(args.log_file)

    run(Path(args.mutation_file),
        Path(args.reference_file),
        Path(args.output_file))
