import os
import pickle
import pandas as pd
from intervaltree import Interval, IntervalTree
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_annotation2(genes_file):
    """Load genes data and rename columns."""
    print("Reading annotation...")
    genes = pd.read_csv(genes_file, sep='\t', header=0, index_col=None)
    genes.rename(columns={'Sequence Id': 'Chromosome'}, inplace=True)
    print(genes)
    genes[['Start', 'Stop']] = genes[['Start', 'Stop']].apply(pd.to_numeric, errors='coerce')
    return genes

def build_and_save_tree(chrom, chrom_df, output_dir):
    """
    Build an IntervalTree for a given chromosome and save to a pickle file.
    Returns (chrom, tree) so we can reconstruct a dictionary if needed.
    """
    itree = IntervalTree(
        Interval(row.Start, row.Stop + 1, row['Locus Tag'])
        for _, row in chrom_df.iterrows()
    )
    file_path = os.path.join(output_dir, f"chromosome_tree_{chrom}.pkl")
    with open(file_path, "wb") as f:
        pickle.dump(itree, f)
    return chrom, itree

def compute_and_save_trees(genes_file, output_dir, threads=4):
    """
    Read genes, group by chromosome, and build/save each IntervalTree in parallel.
    Returns a dict of {chromosome: IntervalTree} for in-memory usage if needed.
    """
    # Read genes
    genes = read_annotation2(genes_file)
    
    # Group the DataFrame by chromosome
    grouped = genes.groupby('Chromosome')
    
    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    results = {}
    # Build & save trees in parallel
    print("Building and saving interval trees (parallelized)...")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # Submit one job per chromosome
        futures = [
            executor.submit(build_and_save_tree, chrom, chrom_df, output_dir)
            for chrom, chrom_df in grouped
        ]
        
        # Collect results as they finish
        for future in as_completed(futures):
            chrom, itree = future.result()
            results[chrom] = itree
    
    print("All interval trees have been built and saved.")
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate or load chromosome trees.")
    parser.add_argument("--genes_file", required=True, help="Gene file to process.")
    parser.add_argument("--output_dir", required=True, help="Directory to save chromosome trees.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use.")
    args = parser.parse_args()
    
    # Parallel execution of building & saving interval trees
    compute_and_save_trees(args.genes_file, args.output_dir, args.threads)
