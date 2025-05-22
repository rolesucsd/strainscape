#!/usr/bin/env python3
"""
Concatenate all per-patient SNV trend/mapped/type tables into one big file.

For each input file it:

1.  extracts the patient ID from the parent directory name
    (…/instrain/<PATIENT_ID>/SNV_filtered_trend_trajectory_mapped_types.txt)
2.  adds a new column **Patient_ID**
3.  appends the rows to the final output

The script streams each file in chunks, so peak RAM stays low even for
multi-gigabyte outputs.

Usage
-----
python combine_trend_mapped_types.py \
       --output combined_trend_mapped_type.txt \
       --inputs file1.tsv file2.tsv … \
       [--chunksize 1_000_000] \
       [--log log.txt]
"""
from __future__ import annotations
import argparse, csv, logging, os, re, sys
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--output", required=True, type=Path,
                   help="Path to combined TSV to write")
    p.add_argument("--inputs", required=True, nargs="+", type=Path,
                   help="List of per-patient TSVs to merge")
    p.add_argument("--chunksize", type=int, default=1_000_000,
                   help="Rows per chunk to read/write")
    p.add_argument("--log", type=Path, help="Optional log file")
    return p.parse_args()

def setup_logging(log_file: Path | None):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(log_file) if log_file else logging.StreamHandler(sys.stderr)
        ],
    )

def extract_patient_id(path: Path) -> str:
    """
    Expect structure …/instrain/<PATIENT_ID>/SNV_filtered_trend_trajectory_mapped_types.txt
    """
    return path.parent.name

def main():
    args = parse_args()
    setup_logging(args.log)

    out_fh = args.output.open("w", buffering=1024 * 1024)
    writer = None

    for tsv in args.inputs:
        pid = extract_patient_id(tsv)
        logging.info(f"→ processing {tsv} (patient {pid})")

        first_chunk = True
        for chunk in pd.read_csv(tsv, sep="\t", chunksize=args.chunksize):
            chunk.insert(0, "Patient_ID", pid)

            # header only on very first write
            header = writer is None and first_chunk
            chunk.to_csv(out_fh, sep="\t", header=header, index=False, mode="a")
            first_chunk = False
            if writer is None:  # pylint: disable=unused-variable
                writer = True  # dummy flag that header has been written
        if first_chunk:  # the file was empty
            logging.warning(f"{tsv} is empty – skipped")

    out_fh.close()
    logging.info(f"✅ combined file written → {args.output}")

if __name__ == "__main__":
    main()
