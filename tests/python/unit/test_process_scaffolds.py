import os
import pandas as pd
import pytest
from snakemake.scripts.process_scaffolds import process_scaffolds

def test_process_scaffolds():
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