#!/usr/bin/env python3
"""
Script to convert a multi-record GenBank file into per-sample FASTQ files
grouped by the 'isolation_source' qualifier of each record.

Adds dummy Phred quality scores to enable importing into QIIME 2 as
SampleData[SequencesWithQuality].

Outputs:
    - Per-sample FASTQ files named after the isolation_source field.
    - A QIIME 2 manifest CSV listing sample-id, path, and direction.

Usage:
    ./genbank_to_fastq_manifest.py -i input.gb -o output_dir -m manifest.csv

Requirements:
    - Biopython (for SeqIO)
"""

import os
import csv
import argparse
from collections import defaultdict

from Bio import SeqIO


def sanitize_string(s):
    """Replace problematic chars with underscores"""
    for ch in [' ', ';', ',', ':']:
        s = s.replace(ch, '_')
    return s


def extract_sample_id(record):
    """
    Extract a sample ID from a GenBank record based on its 'isolation_source' qualifier.

    Args:
        record (SeqRecord): A Biopython SeqRecord from GenBank.

    Returns:
        str: The sample ID with spaces replaced by underscores, or 'UnknownSample'.
    """
    for feat in record.features:
        if "isolation_source" in feat.qualifiers:
            raw_str = feat.qualifiers["isolation_source"][0]
            return sanitize_string(raw_str)
    return "UnknownSample"


def group_and_write_fastq(genbank_file, output_dir, manifest_file):
    """
    Group sequences by sample ID, assign dummy quality scores,
    write per-sample FASTQ files, and generate a QIIME 2 manifest.

    Args:
        genbank_file (str): Path to the input GenBank file.
        output_dir (str): Directory to save FASTQ files.
        manifest_file (str): Path to write QIIME 2 manifest CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    samples = defaultdict(list)

    # Parse GenBank records and group by sample ID
    for record in SeqIO.parse(genbank_file, "genbank"):
        sample_id = extract_sample_id(record)
        # Add dummy quality scores (Phred 40 = ASCII 'I') for compatibility
        record.letter_annotations["phred_quality"] = [40] * len(record.seq)
        samples[sample_id].append(record)

    # Write FASTQ files and manifest
    with open(manifest_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample-id", "absolute-filepath", "direction"])

        for sample_id, records in samples.items():
            fastq_filename = f"{sample_id}.fastq"
            fastq_path = os.path.abspath(os.path.join(output_dir, fastq_filename))

            with open(fastq_path, "w") as out_f:
                SeqIO.write(records, out_f, "fastq")

            writer.writerow([sample_id, fastq_path, "forward"])

    print(f"Wrote {len(samples)} FASTQ files and manifest to: {manifest_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GenBank to per-sample FASTQ files with dummy quality and create QIIME2 manifest."
    )
    parser.add_argument("-i", "--input_gb", required=True, help="Input GenBank file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output FASTQ directory")
    parser.add_argument("-m", "--manifest_file", required=True, help="Output manifest CSV path")

    args = parser.parse_args()
    group_and_write_fastq(args.input_gb, args.output_dir, args.manifest_file)


if __name__ == "__main__":
    main()
