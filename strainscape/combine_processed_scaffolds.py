#!/usr/bin/env python3
"""
Combine processed scaffolds files from multiple patients.

This script concatenates all input TSV files into one large TSV.
It assumes that the Sample column is already present in each input file.

Inputs:
  - A list of processed_scaffolds.txt files
Outputs:
  - A single combined_processed_scaffolds.txt file
"""

import argparse
import pandas as pd
from pathlib import Path

def combine_files(input_files, output_file):
    with open(output_file, 'w') as out_f:
        header_written = False
        for file_path in input_files:
            with open(file_path, 'r') as in_f:
                header = in_f.readline()
                if not header_written:
                    out_f.write(header)
                    header_written = True
                for line in in_f:
                    out_f.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine processed_scaffolds.txt files")
    parser.add_argument('--input-files', nargs='+', required=True,
                        help='Paths to processed_scaffolds.txt files')
    parser.add_argument('--output-file', required=True,
                        help='Path to the combined output file')
    args = parser.parse_args()

    combine_files(args.input_files, args.output_file)
