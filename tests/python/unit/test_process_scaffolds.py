import os
import pandas as pd
import pytest
from strainscape.process_scaffolds import process_scaffolds
import tempfile
from pathlib import Path
import logging
from strainscape.process_scaffolds import process

def test_process_scaffolds():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test files
        scaffold_file = os.path.join(tmpdir, 'scaffold_info.tsv')
        bin_file = os.path.join(tmpdir, 'bin.txt')
        output_file = os.path.join(tmpdir, 'processed_scaffolds.tsv')
        log_file = os.path.join(tmpdir, 'process_scaffolds.log')

        # Create test scaffold data
        pd.DataFrame({
            'scaffold': ['scaffold1', 'scaffold2'],
            'length': [2000, 500],
            'coverage': [10.0, 2.0],
            'breadth': [0.8, 0.3],
            'nucl_diversity': [0.01, 0.02],
            'Sample': ['sample1', 'sample1']
        }).to_csv(scaffold_file, sep='\t', index=False)

        # Create test bin data
        pd.DataFrame({
            'scaffold': ['scaffold1', 'scaffold2'],
            'bin': ['bin1', 'bin2'],
            'Completeness': [90.0, 85.0],
            'Contamination': [5.0, 8.0],
            'Genome_Size': [2000000, 1500000]
        }).to_csv(bin_file, sep='\t', index=False)

        # Run the process_scaffolds function
        process_scaffolds(
            scaffold_file=scaffold_file,
            bin_file=bin_file,
            output_file=output_file,
            min_length=1000,
            min_coverage=5.0,
            min_breadth=0.4,
            min_completeness=50.0,
            max_contamination=10.0,
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

    # Sample input files
    scaffold_file = 'tests/python/unit/data/instrain/M2042/combined/scaffold_info.tsv'
    stb_file = 'tests/python/unit/data/assembly/M2042/assembly.stb'
    metadata_file = 'tests/python/unit/input/metadata/hmp2_metadata_2018-08-20.csv'
    bin_dir = 'tests/python/unit/data/bins/M2042'
    output_file = 'tests/python/unit/data/instrain/M2042/processed_scaffolds_test.tsv'
    log_file = 'tests/python/unit/logs/instrain/M2042/process_scaffolds_test.log'

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