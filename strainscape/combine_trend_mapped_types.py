#!/usr/bin/env python3
"""
Concatenate all per-patient SNV trend/mapped/type tables into one big TSV, preserving all week columns.

For each input file it:
1. extracts the patient ID from the parent directory name
2. adds a new column **Patient_ID**
3. appends the rows to the final output

The script first reads only the headers of each input to compute the union of all columns,
then streams through each file in chunks, reindexes to the full column set, and writes out.

Usage
-----
python combine_trend_mapped_types.py \
       --output combined_trend_mapped_type.txt \
       --inputs file1.tsv file2.tsv … \
       [--chunksize 1000000] \
       [--log log.txt]
"""

import argparse, logging, sys
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--output", required=True, type=Path,
                   help="Path to combined TSV to write")
    p.add_argument("--inputs", required=True, nargs='+', type=Path,
                   help="List of per-patient TSVs to merge")
    p.add_argument("--chunksize", type=int, default=1_000_000,
                   help="Rows per chunk to read/write")
    p.add_argument("--log", type=Path,
                   help="Optional log file (default: stderr)")
    return p.parse_args()


def setup_logging(log_file: Path | None):
    handlers = []
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    else:
        handlers.append(logging.StreamHandler(sys.stderr))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers
    )


def extract_patient_id(path: Path) -> str:
    # expects …/instrain/<PATIENT_ID>/filename
    return path.parent.name


def main():
    args = parse_args()
    setup_logging(args.log)

    # 1) Determine full set of columns
    full_cols = ['Patient_ID']
    for tsv in args.inputs:
        logging.info(f"Reading header from {tsv}")
        header = pd.read_csv(tsv, sep="\t", nrows=0).columns.tolist()
        for col in header:
            if col not in full_cols:
                full_cols.append(col)

    logging.info(f"Final columns: {full_cols}")

    # 2) Open output and write header
    with args.output.open('w', buffering=1024*1024) as out_fh:
        out_fh.write("\t".join(full_cols) + "\n")

        # 3) Stream each file
        for tsv in args.inputs:
            pid = extract_patient_id(tsv)
            logging.info(f"Processing {tsv} (patient {pid})")
            first_chunk = True
            for chunk in pd.read_csv(tsv, sep="\t", chunksize=args.chunksize):
                chunk.insert(0, 'Patient_ID', pid)
                # align to full_cols (missing columns become NaN)
                chunk = chunk.reindex(columns=full_cols)
                # write without header
                chunk.to_csv(out_fh, sep="\t", index=False, header=False)
                first_chunk = False
            if first_chunk:
                logging.warning(f"{tsv} is empty; skipped")

    logging.info(f"✅ combined file written → {args.output}")


if __name__ == '__main__':
    main()
