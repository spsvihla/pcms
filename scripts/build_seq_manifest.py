#!/usr/bin/env python3

import os
import csv
import argparse
from collections import defaultdict

from Bio import SeqIO


def extract_sample_id(record):
    for feat in record.features:
        quals = feat.qualifiers
        if "isolation_source" in quals:
            return quals["isolation_source"][0].replace(" ", "_")
    return "UnknownSample"


def group_and_write_fasta(genbank_file, output_dir, manifest_file):
    os.makedirs(output_dir, exist_ok=True)
    samples = defaultdict(list)

    for record in SeqIO.parse(genbank_file, "genbank"):
        sample_id = extract_sample_id(record)
        samples[sample_id].append(record)

    with open(manifest_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample-id", "absolute-filepath", "direction"])

        for sample_id, records in samples.items():
            fasta_path = os.path.abspath(os.path.join(output_dir, f"{sample_id}.fasta"))
            with open(fasta_path, "w") as out_f:
                SeqIO.write(records, out_f, "fasta")
            writer.writerow([sample_id, fasta_path, "forward"])

    print(f"Wrote {len(samples)} FASTA files and manifest to {manifest_file}")


def main():
    parser = argparse.ArgumentParser(description="Convert GenBank to per-sample FASTA + QIIME2 manifest")

    parser.add_argument(
        '-i', '--input_gb',
        required=True,
        help="GenBank file"
    )
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help="Directory to write FASTA files"
    )
    parser.add_argument(
        '-m', '--manifest_file',
        required=True,
        help="Path to manifest CSV"
    )

    args = parser.parse_args()
    group_and_write_fasta(args.input_gb, args.output_dir, args.manifest_file)


if __name__ == "__main__":
    main()
