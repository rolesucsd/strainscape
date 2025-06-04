import os
import pytest
import pandas as pd
from strainscape.filter_mutations import filter_mutations as filter_mutations_file

def test_filter_mutations():
    # Sample input files
    snv_file = 'tests/python/unit/data/instrain/M2042/combined/snv_info.tsv'
    processed_scaffold_file = 'tests/python/unit/data/instrain/M2042/processed_scaffolds.tsv'
    metadata_file = 'tests/python/unit/input/metadata/hmp2_metadata_2018-08-20.csv'
    output_file = 'tests/python/unit/data/instrain/M2042/filtered_mutations_test.tsv'
    log_file = 'tests/python/unit/logs/instrain/M2042/filter_mutations_test.log'

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Run the filter_mutations function
    filter_mutations_file(
        snv_file=snv_file,
        output_file=output_file,
        metadata_file=metadata_file,
        processed_scaffolds_file=processed_scaffold_file,
        min_coverage=10
    )

    # Verify the output file exists
    assert os.path.exists(output_file), "Output file was not created."

    # Verify the output file contains data
    df = pd.read_csv(output_file, sep='\t')
    assert not df.empty, "Output file is empty."

    # Clean up
    os.remove(output_file) 