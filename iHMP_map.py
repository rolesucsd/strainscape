import os
import pickle
import pandas as pd
from intervaltree import IntervalTree
from joblib import Parallel, delayed
import numpy as np

def read_annotation(genes_file, eggnog_file):
    # Load genes data
    genes = pd.read_csv(genes_file, sep='\t')
    genes.rename(columns={'Sequence Id': 'Chromosome'}, inplace=True)
    genes[['Start', 'Stop']] = genes[['Start', 'Stop']].apply(pd.to_numeric, errors='coerce')
#    genes.dropna(subset=['Start', 'Stop'], inplace=True)

    # Load EggNOG annotations
    eggnog = pd.read_csv(eggnog_file, sep='\t')
    eggnog.rename(columns={'query': 'Locus Tag'}, inplace=True)
    return genes, eggnog

def read_annotation2(genes_file):
    # Load genes data
    genes = pd.read_csv(genes_file, sep='\t')
    genes.rename(columns={'Sequence Id': 'Chromosome'}, inplace=True)
    genes[['Start', 'Stop']] = genes[['Start', 'Stop']].apply(pd.to_numeric, errors='coerce')
#    genes.dropna(subset=['Start', 'Stop'], inplace=True)
    return genes

def load_tree_for_chromosome(chrom, trees_dir):
    file_path = os.path.join(trees_dir, f"chromosome_tree_{chrom}.pkl")
    if os.path.exists(file_path):
        with open(file_path, "rb") as file:
            return pickle.load(file)
    return None

def batch_map_genes(SNV_trajectory, genes, trees_dir):
    results = []
    for chrom, chrom_snvs in SNV_trajectory.groupby('Chromosome'):
        tree = load_tree_for_chromosome(chrom, trees_dir)
        if not tree:
            continue
        # Restrict 'genes' to the chromosome of interest
        chrom_genes = genes[genes['Chromosome'] == chrom]
        
        # Pre-sort by Stop (for finding upstream genes) and by Start (for finding downstream genes)
        chrom_genes_sorted_by_stop = chrom_genes.sort_values('Stop').reset_index(drop=True)
        stops = chrom_genes_sorted_by_stop['Stop'].values
        
        chrom_genes_sorted_by_start = chrom_genes.sort_values('Start').reset_index(drop=True)
        starts = chrom_genes_sorted_by_start['Start'].values

        for _, row in chrom_snvs.iterrows():
            pos = row['Position']
            matches = tree.at(pos)
            
            # --- 1. If there is a genic match in the interval tree
            if matches:
                matched_tag = next(iter(matches)).data
                matched_genes_df = chrom_genes[chrom_genes['Locus Tag'] == matched_tag]

                if not matched_genes_df.empty:
                    matched_gene = matched_genes_df.iloc[0]
                    results.append({
                        'Chromosome': chrom,
                        'Position': pos,
                        'Matched_Locus_Tag': matched_tag,
                        'Matched_Start': matched_gene['Start'],
                        'Matched_Stop': matched_gene['Stop'],
                        'Matched_Strand': matched_gene['Strand'],
                        'Matched_Gene': matched_gene['Gene'],
                        'Matched_Product': matched_gene['Product'],
                        'Matched_DbXrefs': matched_gene['DbXrefs'],
                        'coding': 'genic'
                    })
                else:
                    results.append({
                        'Chromosome': chrom,
                        'Position': pos,
                        'Matched_Locus_Tag': None,
                        'Matched_Start': None,
                        'Matched_Stop': None,
                        'Matched_Strand': None,
                        'Matched_Gene': None,
                        'Matched_Product': None,
                        'Matched_DbXrefs': None,
                        'coding': 'genic'
                    })
            # --- 2. Intergenic case
            else:
                result = {
                    'Chromosome': chrom,
                    'Position': pos,
                    'coding': 'intergenic'
                }
                
                # Find the index of the first gene in 'stops' whose 'Stop' >= pos
                pos_idx = np.searchsorted(stops, pos)
                # The upstream gene would be the one immediately before this index
                upstream_idx = pos_idx - 1
                
                if 0 <= upstream_idx < len(chrom_genes_sorted_by_stop):
                    upstream_gene = chrom_genes_sorted_by_stop.iloc[[upstream_idx]]
                else:
                    # define an empty DataFrame if no upstream gene
                    upstream_gene = pd.DataFrame()
                
                # Similarly, find the downstream gene using the 'starts' array
                # The first gene whose 'Start' > pos is at index = searchsorted(starts, pos, 'right')
                # but typically we use 'searchsorted(starts, pos, side="left")' to find the insertion point
                down_idx = np.searchsorted(starts, pos, side="left")
                if 0 <= down_idx < len(chrom_genes_sorted_by_start):
                    downstream_gene = chrom_genes_sorted_by_start.iloc[[down_idx]]
                else:
                    downstream_gene = pd.DataFrame()

                # If upstream_gene is not empty, update the result dict
                if not upstream_gene.empty:
                    result.update({
                        'upstream_Locus_Tag': upstream_gene['Locus Tag'].values[0],
                        'upstream_Start': upstream_gene['Start'].values[0],
                        'upstream_Stop': upstream_gene['Stop'].values[0],
                        'upstream_Gene': upstream_gene['Gene'].values[0],
                        'upstream_Product': upstream_gene['Product'].values[0],
                        'upstream_DbXrefs': upstream_gene['DbXrefs'].values[0]
                    })

                # If downstream_gene is not empty, update the result dict
                if not downstream_gene.empty:
                    result.update({
                        'downstream_Locus_Tag': downstream_gene['Locus Tag'].values[0],
                        'downstream_Start': downstream_gene['Start'].values[0],
                        'downstream_Stop': downstream_gene['Stop'].values[0],
                        'downstream_Gene': downstream_gene['Gene'].values[0],
                        'downstream_Product': downstream_gene['Product'].values[0],
                        'downstream_DbXrefs': downstream_gene['DbXrefs'].values[0]
                    })

                results.append(result)

    return pd.DataFrame(results)
