import pandas as pd
import tempfile
import os
from pathlib import Path
from strainscape.map_genes import map_genes_fast

def test_map_genes_fast():
    # Create dummy mutation, trend, and gene files
    with tempfile.TemporaryDirectory() as tmpdir:
        trend_file = Path(tmpdir) / "trend.tsv"
        gene_file = Path(tmpdir) / "genes.tsv"
        out_file = Path(tmpdir) / "out.tsv"

        # Write minimal trend file
        pd.DataFrame({
            'scaffold': ['chr1'],
            'position': [100],
            'ref_base': ['A'],
            'new_base': ['G']
        }).to_csv(trend_file, sep='\t', index=False)

        # Write minimal gene file with Bakta header
        with open(gene_file, 'w') as f:
            f.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n')
            f.write('chr1\tCDS\t50\t150\t+\tLT1\tGENE1\tproduct\tNA\n')

        # Run the function
        map_genes_fast(trend_file, gene_file, out_file)

        # Check output file exists and is not empty
        assert out_file.exists()
        df = pd.read_csv(out_file, sep='\t')
        assert not df.empty 