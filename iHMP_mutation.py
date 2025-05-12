import pandas as pd
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_genic_vs_intergenic_and_mutations(df, coding_col='coding', mutation_col='Mutation_Type', save_fig=None):
    """
    Plots two subplots:
      1) Distribution of 'genic' vs 'intergenic' 
      2) Distribution of mutation types (e.g. Silent, Missense, Nonsense)

    :param df: pd.DataFrame containing at least the columns [coding_col, mutation_col].
    :param coding_col: Name of the column that labels a SNP as genic vs intergenic.
    :param mutation_col: Name of the column with mutation types (e.g. Silent, Missense).
    :param save_fig: str or None. If provided, saves the figure to this path.
    """
    # Prepare data for genic vs intergenic
    coding_counts = df[coding_col].value_counts(dropna=False)

    # Prepare data for mutation types
    mutation_counts = df[mutation_col].value_counts(dropna=False)

    # Create figure and axes
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    ax1, ax2 = axes

    # 1) Bar chart (or pie) of genic vs intergenic
    sns.barplot(
        x=coding_counts.index.astype(str),
        y=coding_counts.values,
        ax=ax1,
        palette='Set2'
    )
    ax1.set_title("Genic vs. Intergenic")
    ax1.set_xlabel("Category")
    ax1.set_ylabel("Count of SNPs")

    # If you prefer a pie chart, uncomment these lines instead and comment out barplot:
    # ax1.pie(coding_counts.values, labels=coding_counts.index, autopct='%1.1f%%', startangle=140)
    # ax1.set_title("Genic vs. Intergenic")

    # 2) Bar chart (or pie) of mutation types
    sns.barplot(
        x=mutation_counts.index.astype(str),
        y=mutation_counts.values,
        ax=ax2,
        palette='Set3'
    )
    ax2.set_title("Mutation Types")
    ax2.set_xlabel("Type")
    ax2.set_ylabel("Count of SNPs")

    # If you prefer a pie chart, uncomment:
    # ax2.pie(mutation_counts.values, labels=mutation_counts.index, autopct='%1.1f%%', startangle=140)
    # ax2.set_title("Mutation Types")

    plt.tight_layout()

    if save_fig:
        plt.savefig(save_fig+"/mutation.png", dpi=600)
    else:
        plt.show()


def get_mutation_type(row, sequences):
    """
    Determine the mutation type (Silent, Missense, Nonsense) for a single SNV row.

    Assumptions:
    1. The 'old_base' (ref_base) and 'new_base' are in the same (forward) orientation
       as 'sequences[chrom]', i.e., the plus-strand reference.
    2. gene_start and gene_stop are 1-based coordinates on the plus strand of the chromosome.
       Even if the gene is on the minus strand, gene_start < gene_stop (the typical GFF convention).
    3. position is also 1-based on the plus strand.

    If your data differs from these assumptions, adjust accordingly.
    """
    # If it's not genic or doesn't have a valid Matched_Start/Stop, skip
    if row['coding'] != 'genic' or pd.isna(row['Matched_Start']) or pd.isna(row['Matched_Stop']):
        return pd.Series({'Mutation_Type': None})

    # Extract relevant columns
    chrom       = row['Chromosome']                # Chrom/scaffold name
    position    = int(row['Position'])             # 1-based coordinate
    gene_start  = int(row['Matched_Start'])        # 1-based
    gene_stop   = int(row['Matched_Stop'])         # 1-based
    strand      = row['Matched_Strand']            # '+' or '-'
    old_base    = row['ref_base'].upper()
    new_base    = row['new_base'].upper()

    # Get the chromosome sequence in forward orientation
    chrom_seq = sequences.get(chrom)
    if chrom_seq is None:
        return pd.Series({'Mutation_Type': 'Gene sequence not found'})

    gene_len = gene_stop - gene_start + 1
    if gene_len <= 0:
        return pd.Series({'Mutation_Type': 'Invalid gene range'})

    # Extract the gene region in the *forward orientation* of the chromosome
    #    Convert 1-based -> 0-based for slicing
    gene_substr = chrom_seq[gene_start - 1 : gene_stop]  # forward substring

    # Compute the offset of 'position' within that substring
    pos_index = position - gene_start + 1

    # Check bounds
    if pos_index < 0 or pos_index >= gene_len:
        return pd.Series({'Mutation_Type': f'Position {position} out of gene bounds'})

    # Current_base in gene_substr (forward orientation)
    current_base = gene_substr[pos_index].upper()

    
    # Compare the forward-orientation base to old_base
    #    If your 'old_base' is truly forward-based, then they must match directly
    #    for both + and - strand genes (since gene_substr is forward).
    
    if current_base != old_base:
        return pd.Series({'Mutation_Type': f'Reference base mismatch (expected {old_base}, got {current_base})'})

    # Construct the mutated substring in the forward orientation
    mutated_substr_list = list(gene_substr)
    mutated_substr_list[pos_index] = new_base
    mutated_substr_fwd = "".join(mutated_substr_list)

    # Identify codon boundaries
    codon_start = (pos_index // 3) * 3
    original_codon_fwd = gene_substr[codon_start : codon_start + 3]
    mutated_codon_fwd  = mutated_substr_fwd[codon_start : codon_start + 3]

    # If the gene is minus strand, the *actual transcribed codon* is the reverse complement
    if strand == '-':
        original_codon = str(Seq(original_codon_fwd).reverse_complement())
        mutated_codon  = str(Seq(mutated_codon_fwd).reverse_complement())
    else:
        original_codon = original_codon_fwd
        mutated_codon  = mutated_codon_fwd

    # Translate
    original_aa = str(Seq(original_codon).translate())
    mutated_aa  = str(Seq(mutated_codon).translate())

    if original_aa == mutated_aa:
        return pd.Series({'Mutation_Type': 'Silent'})
    elif mutated_aa == '*' and original_aa != '*':
        return pd.Series({'Mutation_Type': 'Nonsense'})
    else:
        return pd.Series({'Mutation_Type': 'Missense'})
