#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ===============================
# Config: Set variables here
# ===============================

DATA="${DATA:-/path/to/data}"

ACCESSIONS="JN427016:JN539989"
PROJECT="guerrero_negro"

DATA_DIR="$DATA/$PROJECT"
BUILD_DIR="$DATA_DIR/build"

REF_SEQ_FA="$DATA/greengenes/gg_13_8_otus/rep_set/97_otus.fasta"
REF_CHIMERA_FA="$DATA/gold_nodup.fa"

POSTFIX="${ACCESSIONS//:/_}"

# Files derived from accessions
ACC_LIST="$DATA_DIR/accessions_${POSTFIX}.txt"
GENBANK="$DATA_DIR/query_seqs_${POSTFIX}.gb"
FASTQ_DIR="$BUILD_DIR/query_seqs_${POSTFIX}"
MANIFEST="$BUILD_DIR/query_manifest_${POSTFIX}.csv"
SEQS_QZA="$BUILD_DIR/query_seqs_${POSTFIX}.qza"
TABLE_QZA="$BUILD_DIR/query_table_${POSTFIX}.qza"
DEREP_SEQS_QZA="$BUILD_DIR/query_seqs_derep_${POSTFIX}.qza"

CLUST_PREFIX="gg_13_8_97ref_97clust"
CLUST_DIR="$BUILD_DIR/$CLUST_PREFIX"
CLUST_TABLE_QZA="$CLUST_DIR/${CLUST_PREFIX}_table.qza"
CLUST_TABLE_TSV="$CLUST_DIR/${CLUST_PREFIX}_table.tsv"
CLUST_SEQS_QZA="$CLUST_DIR/${CLUST_PREFIX}_seqs.qza"
CLUST_SEQS_FA="$CLUST_DIR/${CLUST_PREFIX}_seqs.fa"
NEW_REFS_QZA="$CLUST_DIR/${CLUST_PREFIX}_new_refs.qza"
NEW_REFS_FA="$CLUST_DIR/${CLUST_PREFIX}_new_refs.fa"

# Chimera filtering outputs
FILTERED_PREFIX="${CLUST_PREFIX}_filtered"
FILTERED_DIR="$BUILD_DIR/${FILTERED_PREFIX}"
CHIM_TABLE_QZA="$FILTERED_DIR/${FILTERED_PREFIX}_table.qza"
CHIM_TABLE_TSV="$DATA_DIR/${FILTERED_PREFIX}_table.tsv"
CHIM_SEQS_QZA="$FILTERED_DIR/${FILTERED_PREFIX}_seqs.qza"
CHIM_SEQS_FA="$DATA_DIR/${FILTERED_PREFIX}_seqs.fa"
# CHIM_IDS="$FILTERED_DIR/${CLUST_PREFIX}_chimera_ids.txt"

# ===============================
# Helper: Print step header
# ===============================

step() {
  echo -e "\n=== Step $1: $2 ==="
}

# ===============================
# Main workflow
# ===============================

cd ./qiime2-utils

mkdir -p "$BUILD_DIR" "$FASTQ_DIR" "$CLUST_DIR" "$FILTERED_DIR"

step 1 "Build accession list"
./fetch_accession_list.sh \
  --i-accession "$ACCESSIONS" \
  --o-outdir "$DATA_DIR"

step 2 "Fetch GenBank files"
./fetch_ncbi_genbank_files.sh \
  --i-accession-list "$ACC_LIST" \
  --o-genbank "$GENBANK" \
  --batch-size 100

step 3 "Convert GenBank to FASTQ and build manifest"
python3 build_fastq_manifest.py \
  --i-genbank "$GENBANK" \
  --o-outdir "$FASTQ_DIR" \
  --o-manifest "$MANIFEST"

step 4 "Import query sequences into QIIME2"
./import_query_seqs_from_manifest.sh \
  --i-manifest "$MANIFEST" \
  --o-sequences "$SEQS_QZA"

step 5 "Dereplicate sequences"
./dereplicate_seqs.sh \
  --i-sequences "$SEQS_QZA" \
  --o-table "$TABLE_QZA" \
  --o-sequences "$DEREP_SEQS_QZA"

step 6 "Open-reference clustering"
./open_ref_cluster_seqs.sh \
  --i-ref-fasta "$REF_SEQ_FA" \
  --i-table "$TABLE_QZA" \
  --i-query-qza "$DEREP_SEQS_QZA" \
  --i-identity 0.97 \
  --prefix "$CLUST_PREFIX" \
  --o-outdir "$CLUST_DIR"

step 7 "Convert table to TSV"
./export_table_to_tsv_from_qza.sh \
  -q "$CLUST_TABLE_QZA" \
  -t "$CLUST_TABLE_TSV"

step 8 "Export new reference sequences to FASTA"
./export_seqs_to_fasta_from_qza.sh \
  -q "$CLUST_SEQS_QZA" \
  -f "$CLUST_SEQS_FA"

step 9 "Export de-novo clustered sequences to FASTA"
./export_seqs_to_fasta_from_qza.sh \
  -q "$NEW_REFS_QZA" \
  -f "$NEW_REFS_FA"

step 10 "Filter chimeric sequences"
./filter_chimeras.sh \
  --i-sequences "$CLUST_SEQS_QZA" \
  --i-table "$CLUST_TABLE_QZA" \
  --i-ref-chimera-fasta "$REF_CHIMERA_FA" \
  --prefix "$FILTERED_PREFIX" \
  --o-outdir "$FILTERED_DIR"

step 11 "Convert chimera-filtered table to TSV"
./export_table_to_tsv_from_qza.sh \
  -q "$CHIM_TABLE_QZA" \
  -t "$CHIM_TABLE_TSV"

step 12 "Export chimera-filtered new reference sequences to FASTA"
./export_seqs_to_fasta_from_qza.sh \
  -q "$CHIM_SEQS_QZA" \
  -f "$CHIM_SEQS_FA"

echo -e "\nPipeline complete!"
echo "Chimera-filtered OTU Table: $CHIM_TABLE_TSV"
echo "Chimera-filtered OTU Sequences: $CHIM_SEQS_FA"
echo "New OTUs Added: $NEW_REFS_FA"
