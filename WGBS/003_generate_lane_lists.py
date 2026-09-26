#!/usr/bin/env python3
"""
01_generate_lane_lists.py

Parse an EGA sample CSV metafile and produce two output files:

  1. all_samples_lanes.txt
     All samples and their lane prefixes in one file.
     Format:
         sample_alias  lane_prefix
         tumor_4121361 WGBS-PAIRED-tumor_4121361-SN509_0187_BC130GACXX_5
         tumor_4121361 WGBS-PAIRED-tumor_4121361-SN1007_0108_AC0T0UACXX_1
         ...

  2. completed_samples.txt
     One line per sample that has ALL expected lane BAMs already present
     under the BWA directory (*_filtered.bam).
     Format:
         sample_alias  n_expected  n_found  bam_path1,bam_path2,...

Usage:
    python3 01_generate_lane_lists.py \
        --csv   sample_file_with_checksums_Apr23_2026.csv \
        --bwa   /../BWA \
        --outdir .
"""

import csv
import argparse
import os
import glob
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Generate lane list and completed-sample list from EGA CSV metadata."
    )
    parser.add_argument("--csv",     required=True,
                        help="Input CSV metafile (EGA sample/file manifest)")
    parser.add_argument("--bwa",     required=True,
                        help="Path to BWA output directory containing *_filtered.bam files")
    parser.add_argument("--outdir",  default=".",
                        help="Directory for output txt files (default: current directory)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Parse CSV: collect R1 lane prefixes per sample_alias
    # ------------------------------------------------------------------
    samples = defaultdict(list)
    with open(args.csv, newline='') as fh:
        for row in csv.DictReader(fh):
            fname = row['file_name'].strip()
            alias = row['sample_alias'].strip()
            if fname.endswith('-R1.fastq.gz'):
                prefix = fname.replace('-R1.fastq.gz', '')
                samples[alias].append(prefix)

    # ------------------------------------------------------------------
    # 2. Write all_samples_lanes.txt  (all samples, one lane per line)
    # ------------------------------------------------------------------
    lanes_path = os.path.join(args.outdir, "all_samples_lanes.txt")
    with open(lanes_path, 'w') as fh:
        fh.write("sample_alias\tlane_prefix\n")
        for alias in sorted(samples):
            for prefix in samples[alias]:
                fh.write(f"{alias}\t{prefix}\n")

    total_lanes = sum(len(v) for v in samples.values())
    print(f"[1] Lane list written : {lanes_path}")
    print(f"    {len(samples)} samples, {total_lanes} lanes total")

    # ------------------------------------------------------------------
    # 3. Scan BWA directory for *_filtered.bam and check completion
    # ------------------------------------------------------------------
    # Index all found filtered BAMs by sample alias
    found_bams = defaultdict(list)
    for bam in sorted(glob.glob(os.path.join(args.bwa, "*_filtered.bam"))):
        basename = os.path.basename(bam)   # e.g. WGBS-PAIRED-tumor_4121361-SN509_..._filtered.bam
        # Match against known sample aliases (longest match wins)
        matched_alias = None
        for alias in samples:
            tag = f"WGBS-PAIRED-{alias}-"
            if basename.startswith(tag):
                matched_alias = alias
                break
        if matched_alias:
            found_bams[matched_alias].append(bam)

    # ------------------------------------------------------------------
    # 4. Write completed_samples.txt
    #    A sample is "completed" when n_found == n_expected
    # ------------------------------------------------------------------
    completed_path = os.path.join(args.outdir, "completed_samples.txt")
    completed_count = 0
    incomplete_count = 0

    with open(completed_path, 'w') as fh:
        fh.write("sample_alias\tn_expected\tn_found\tstatus\tbam_files\n")
        for alias in sorted(samples):
            n_expected = len(samples[alias])
            found      = sorted(found_bams.get(alias, []))
            n_found    = len(found)
            status     = "COMPLETE" if n_found == n_expected else "INCOMPLETE"
            bam_list   = ",".join(found) if found else "NA"
            fh.write(f"{alias}\t{n_expected}\t{n_found}\t{status}\t{bam_list}\n")

            if status == "COMPLETE":
                completed_count += 1
            else:
                incomplete_count += 1
                # Print which lanes are missing
                found_prefixes = {
                    os.path.basename(b).replace('_filtered.bam', '') for b in found
                }
                missing = [p for p in samples[alias] if p not in found_prefixes]
                print(f"  [INCOMPLETE] {alias}: {n_found}/{n_expected} BAMs found")
                for m in missing:
                    print(f"    missing: {m}_filtered.bam")

    print(f"\n[2] Completed sample list written : {completed_path}")
    print(f"    COMPLETE   : {completed_count} samples")
    print(f"    INCOMPLETE : {incomplete_count} samples")


if __name__ == "__main__":
    main()
