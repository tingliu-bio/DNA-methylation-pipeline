#!/usr/bin/env python3
"""
make_samples_csv.py
───────────────────
Convert the existing sample manifest (tab-separated: Patient_ID, Time_Point,
hpc_path_R1) into a Snakemake-compatible samples.csv with one row per lane.

Parses flowcell ID and lane number directly from the R1 filename.
R2 path is derived by replacing '_R1.fastq.gz' with '_R2.fastq.gz'.

Input format (tab-separated, with header):
    Patient_ID  Time_Point  hpc_path_R1  [source_batch]
    SAMPLE_A_T1 T1          /cluster/.../..._A225FC7LT3_1_CTGTACCA-AGTCTGTG_R1.fastq.gz

Output: samples.csv
    sample, lane, flowcell, R1, R2

Usage:
    python make_samples_csv.py --input manifest.txt --output samples.csv
"""

import argparse
import os
import re
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert sample manifest to Snakemake samples.csv"
    )
    p.add_argument("--input",  required=True,
                   help="Input manifest file (tab-separated, with header)")
    p.add_argument("--output", default="samples.csv",
                   help="Output CSV file (default: samples.csv)")
    return p.parse_args()


def parse_flowcell_lane(r1_path):
    """
    Extract flowcell ID and lane number from an R1 filename.

    Expected filename pattern (anywhere in the name):
        ..._<FLOWCELL>_<LANE>_<BARCODE>_R1.fastq.gz
    where:
        FLOWCELL = alphanumeric string (e.g. A225FC7LT3)
        LANE     = 1-2 digit integer (e.g. 1, 2, ..., 8)
        BARCODE  = dual-index string (e.g. CTGTACCA-AGTCTGTG)

    Returns:
        (flowcell, lane_str)  e.g. ('A225FC7LT3', 'L001')
    Raises:
        ValueError if the pattern is not found.
    """
    fname = os.path.basename(r1_path)

    # Match: _FLOWCELL_LANE_BARCODE_R1.fastq.gz at end of filename
    # FLOWCELL: uppercase letters and digits
    # LANE: 1-2 digits
    # BARCODE: uppercase letters, digits and hyphens
    pattern = r'_([A-Z0-9]+)_(\d{1,2})_[A-Z]+-[A-Z]+_R1\.fastq\.gz$'
    m = re.search(pattern, fname)

    if not m:
        raise ValueError(
            f"Cannot parse flowcell/lane from filename: {fname}\n"
            f"Expected pattern: _FLOWCELL_LANE_BARCODE_R1.fastq.gz"
        )

    flowcell = m.group(1)
    lane     = f"L{m.group(2).zfill(3)}"   # e.g. '1' -> 'L001'
    return flowcell, lane


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input file not found: {args.input}")

    rows = []
    errors = []

    with open(args.input) as fh:
        header = fh.readline().rstrip("\n\r")   # skip header line
        expected_cols = ["Patient_ID", "Time_Point", "hpc_path_R1"]
        # basic header check
        if not all(c in header for c in expected_cols):
            print(f"WARNING: unexpected header: {header}", file=sys.stderr)

        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n\r")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                print(f"WARNING: line {lineno} has fewer than 3 columns, skipping",
                      file=sys.stderr)
                continue

            patient_id = parts[0].strip()
            time_point = parts[1].strip()
            r1_path    = parts[2].strip()

            # sample name combines patient ID and time point
            sample = f"{patient_id}_{time_point}"

            # derive R2 path from R1
            if "_R1.fastq.gz" not in r1_path:
                errors.append(f"Line {lineno}: R1 path does not contain '_R1.fastq.gz': {r1_path}")
                continue
            r2_path = r1_path.replace("_R1.fastq.gz", "_R2.fastq.gz")

            # parse flowcell and lane from filename
            try:
                flowcell, lane = parse_flowcell_lane(r1_path)
            except ValueError as e:
                errors.append(f"Line {lineno}: {e}")
                continue

            rows.append({
                "sample":   sample,
                "lane":     lane,
                "flowcell": flowcell,
                "R1":       r1_path,
                "R2":       r2_path,
            })

    # report parsing errors before writing output
    if errors:
        print("\nParsing errors:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        if not rows:
            sys.exit("ERROR: no rows successfully parsed -- check input format")

    df = pd.DataFrame(rows, columns=["sample", "lane", "flowcell", "R1", "R2"])

    # summary
    print(f"Parsed {len(df)} units across {df['sample'].nunique()} samples")
    print(f"\nLanes per sample:")
    print(df.groupby("sample")["lane"].count().to_string())

    # sanity check: warn if R2 files don't exist
    missing_r2 = [r for r in df["R2"] if not os.path.exists(r)]
    if missing_r2:
        print(f"\nWARNING: {len(missing_r2)} R2 files not found on disk:")
        for f in missing_r2[:5]:
            print(f"  {f}")
        if len(missing_r2) > 5:
            print(f"  ... and {len(missing_r2) - 5} more")

    df.to_csv(args.output, index=False)
    print(f"\nOutput written to: {args.output}")


if __name__ == "__main__":
    main()
