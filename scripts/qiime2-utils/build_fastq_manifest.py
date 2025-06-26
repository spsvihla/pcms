#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from collections import defaultdict

from Bio import SeqIO


# ANSI color codes
CYAN = '\033[0;36m'
GREEN = '\033[0;32m'
NC = '\033[0m'  # No color


def sanitize_string(s):
    """Replace problematic chars with underscores"""
    for ch in [' ', ';', ',', ':']:
        s = s.replace(ch, '_')
    return s


def extract_sample_id(record):
    """Extract a sample ID from a GenBank record based on its 'isolation_source' qualifier."""
    for feat in record.features:
        if "isolation_source" in feat.qualifiers:
            raw_str = feat.qualifiers["isolation_source"][0]
            return sanitize_string(raw_str)
    return "UnknownSample"


def group_and_write_fastq(genbank_file, output_dir, manifest_file):
    """
    Group sequences by sample ID, assign dummy quality scores,
    write per-sample FASTQ files, and generate a QIIME 2 manifest.
    """
    os.makedirs(output_dir, exist_ok=True)
    samples = defaultdict(list)

    print(f"{CYAN}Parsing GenBank file: {genbank_file}{NC}")

    total_records = 0
    for record in SeqIO.parse(genbank_file, "genbank"):
        total_records += 1
        sample_id = extract_sample_id(record)
        # Add dummy quality scores (Phred 40 = ASCII 'I') for compatibility
        record.letter_annotations["phred_quality"] = [40] * len(record.seq)
        samples[sample_id].append(record)

        if total_records % 100 == 0:
            sys.stdout.write(f"\r{CYAN}Processed {total_records} records...{NC}")

    print(f"{CYAN}Total records processed: {total_records}{NC}")
    print(f"{CYAN}Writing FASTQ files to directory: {output_dir}{NC}")

    with open(manifest_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample-id", "absolute-filepath", "direction"])

        for sample_id, records in samples.items():
            fastq_filename = f"{sample_id}.fastq"
            fastq_path = os.path.abspath(os.path.join(output_dir, fastq_filename))

            with open(fastq_path, "w") as out_f:
                SeqIO.write(records, out_f, "fastq")

            writer.writerow([sample_id, fastq_path, "forward"])

    print(f"{GREEN}Finished writing {len(samples)} FASTQ files and manifest:{NC} {manifest_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GenBank to per-sample FASTQ files with dummy quality and create QIIME2 manifest."
    )
    parser.add_argument("--i-genbank", required=True, help="Input GenBank file")
    parser.add_argument("--o-outdir", required=True, help="Output FASTQ directory")
    parser.add_argument("--o-manifest", required=True, help="Output manifest CSV path")

    args = parser.parse_args()
    group_and_write_fastq(args.i_genbank, args.o_outdir, args.o_manifest)


if __name__ == "__main__":
    main()
