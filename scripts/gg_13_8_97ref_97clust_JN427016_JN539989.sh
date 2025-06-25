#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ===============================
# Configuration: Set variables here
# ===============================

# Base data directory (adjust or pass via env)
DATA="${DATA:-/path/to/data}"

# Define dataset parameters
ACCESSION_RANGE="JN427016:JN539989"
DATASET_NAME="guerrero_negro"
BUILD_SUBDIR="build"

# Paths
DATASET_DIR="$DATA/$DATASET_NAME"
BUILD_DIR="$DATASET_DIR/$BUILD_SUBDIR"
REF_SEQS="$DATA/greengenes/gg_13_8_otus/rep_set/97_otus.fasta"

# Derived filenames
ACCESSION_LIST="$DATASET_DIR/accessions_${ACCESSION_RANGE//:/_}.txt"
GENBANK_FILE="$DATASET_DIR/query_seqs_${ACCESSION_RANGE//:/_}.gb"
QUERY_FASTQ_DIR="$BUILD_DIR/query_seqs_${ACCESSION_RANGE//:/_}_fastq"
QUERY_MANIFEST="$BUILD_DIR/query_manifest_${ACCESSION_RANGE//:/_}.csv"
QUERY_SEQS_QZA="$BUILD_DIR/query_seqs_${ACCESSION_RANGE//:/_}.qza"
QUERY_TABLE_QZA="$BUILD_DIR/query_table_${ACCESSION_RANGE//:/_}.qza"
QUERY_DEREP_SEQS_QZA="$BUILD_DIR/query_seqs_derep_${ACCESSION_RANGE//:/_}.qza"
CLUSTER_PREFIX="gg_13_8_97ref_97clust_${ACCESSION_RANGE//:/_}"
CLUSTER_TABLE_QZA="$BUILD_DIR/${CLUSTER_PREFIX}_table.qza"
CLUSTER_TABLE_TSV="$DATASET_DIR/${CLUSTER_PREFIX}_table.tsv"

# ===============================
# Helper: Print step header
# ===============================
step() {
  echo -e "\n=== Step $1: $2 ==="
}

# ===============================
# Main workflow
# ===============================

cd ./open_ref_cluster_from_genbank/

mkdir -p "$BUILD_DIR"
mkdir -p "$QUERY_FASTQ_DIR"

step 1 "Build accession list"
./build_accession_list.sh -a "$ACCESSION_RANGE" -o "$DATASET_DIR"

step 2 "Fetch GenBank files"
./fetch_ncbi_genbank_files.sh -i "$ACCESSION_LIST" -o "$GENBANK_FILE"

step 3 "Convert GenBank to FASTQ and build manifest"
python3 build_seq_manifest.py -i "$GENBANK_FILE" -o "$QUERY_FASTQ_DIR" -m "$QUERY_MANIFEST"

step 4 "Import query sequences into QIIME2"
./import_query_seqs_from_manifest.sh -m "$QUERY_MANIFEST" -o "$QUERY_SEQS_QZA"

step 5 "Dereplicate sequences"
./dereplicate_seqs.sh -i "$QUERY_SEQS_QZA" -t "$QUERY_TABLE_QZA" -s "$QUERY_DEREP_SEQS_QZA"

step 6 "Open-reference clustering"
./open_ref_cluster.sh -r "$REF_SEQS" -t "$QUERY_TABLE_QZA" -q "$QUERY_DEREP_SEQS_QZA" -i 0.97 -p "$CLUSTER_PREFIX" -o "$BUILD_DIR"

step 7 "Convert QIIME2 artifact to TSV"
./convert_qza_table2csv.sh -i "$CLUSTER_TABLE_QZA" -o "$CLUSTER_TABLE_TSV"

echo -e "\nPipeline complete! Results in $CLUSTER_TABLE_TSV"
