import os
import pandas as pd
import pytest
from snakemake.scripts.process_scaffolds import process_scaffolds
import tempfile
from pathlib import Path
import logging
from snakemake.scripts.process_scaffolds import process

def test_process_scaffolds():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        # Create dummy bin_dir with a .fa file
        bin_dir = tmpdir / "bins"
        bin_dir.mkdir()
        fa_file = bin_dir / "bin1.fa"
        with open(fa_file, 'w') as f:
            f.write(">scaf1\nATGC\n")
        # Create dummy scaffold_info.tsv
        scaff_file = tmpdir / "scaffold_info.tsv"
        pd.DataFrame({
            'scaffold': ['scaf1'],
            'length': [2000],
            'coverage': [10.0],
            'breadth': [0.8]
        }).to_csv(scaff_file, sep='\t', index=False)
        # Create dummy stb.tsv
        stb_file = tmpdir / "assembly.stb"
        pd.DataFrame({
            'scaffold': ['scaf1'],
            'Sample': ['S1']
        }).to_csv(stb_file, sep='\t', index=False, header=False)
        # Create dummy metadata.csv
        meta_file = tmpdir / "meta.csv"
        pd.DataFrame({
            'External.ID': ['S1'],
            'week_num': [1],
            'Participant ID': ['P1'],
            'sex': ['F'],
            'diagnosis': ['CD'],
            'Height': [160],
            'Weight': [60],
            'BMI': [23],
            'fecalcal_ng_ml': [100],
            'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)': ['No'],
            'Antibiotics': ['No'],
            'Immunosuppressants (e.g. oral corticosteroids)': ['No']
        }, columns=[
            'External.ID',
            'week_num',
            'Participant ID',
            'sex',
            'diagnosis',
            'Height',
            'Weight',
            'BMI',
            'fecalcal_ng_ml',
            'Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)',
            'Antibiotics',
            'Immunosuppressants (e.g. oral corticosteroids)'
        ]).to_csv(meta_file, index=False)
        # Output file
        out_file = tmpdir / "out.tsv"
        # Run the function
        process(scaff_file, stb_file, meta_file, bin_dir, out_file, log=logging.getLogger("test"))
        # Check output file exists and is not empty
        assert out_file.exists()
        df = pd.read_csv(out_file, sep='\t')
        assert not df.empty

    # Sample input files
    scaffold_file = 'data/instrain/M2042/combined/scaffold_info.tsv'
    stb_file = 'data/assembly/M2042/assembly.stb'
    metadata_file = 'input/metadata/hmp2_metadata_2018-08-20.csv'
    bin_dir = 'data/bins/M2042'
    output_file = 'data/instrain/M2042/processed_scaffolds_test.tsv'
    log_file = 'logs/instrain/M2042/process_scaffolds_test.log'

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Run the process_scaffolds function
    process_scaffolds(
        scaffold_file=scaffold_file,
        stb_file=stb_file,
        metadata_file=metadata_file,
        bin_dir=bin_dir,
        output_file=output_file,
        min_length=1000,
        min_coverage=5.0,
        min_breadth=0.4,
        threads=4,
        chunksize=10000,
        log_file=log_file
    )

    # Verify the output file exists
    assert os.path.exists(output_file), "Output file was not created."

    # Verify the output file contains data
    df = pd.read_csv(output_file, sep='\t')
    assert not df.empty, "Output file is empty."

    # Clean up
    os.remove(output_file)
    os.remove(log_file) 