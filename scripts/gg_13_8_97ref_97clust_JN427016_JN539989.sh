#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Config ---

ACCESSIONS="JN427016:JN539989"
PROJECT="guerrero_negro"

DATA_DIR="$DATA/$PROJECT"
BUILD_DIR="$DATA_DIR/build"

# Reference tree, sequences, and alignment
THRESH=0.05
SAFE_THRESH="${THRESH//./_}"
REF_TREE_NWK="$DATA/greengenes/gg_13_8_otus/trees/97_otus_unannotated.nwk"
REF_TREE_PRUNED_NWK="$DATA/greengenes/gg_13_8_otus/trees/97_otus_unannotated_${SAFE_THRESH}_pruned.nwk"
REF_SEQS_FA="$DATA/greengenes/gg_13_8_otus/rep_set/97_otus.fasta"
REF_SEQS_REDUCED_FA="$DATA/greengenes/gg_13_8_otus/rep_set/97_otus_${SAFE_THRESH}_reduced.fasta"
REF_SEQS_ALIGNED_FA="$DATA/greengenes/gg_13_8_otus/rep_set_aligned/97_otus.fasta"
REF_SEQS_ALIGNED_REDUCED_FA="$DATA/greengenes/gg_13_8_otus/rep_set_aligned/97_otus_${SAFE_THRESH}_reduced.fasta"
REF_CHIMERA_FA="$DATA/gold_nodup.fa"

TRIMMED_REF_SEQ_IDS_TXT="${REF_SEQS_ALIGNED_FA%.fasta}"
TRIMMED_REF_SEQ_IDS_TXT="${TRIMMED_REF_SEQ_IDS_TXT%.fa}"
TRIMMED_REF_SEQ_IDS_TXT="${TRIMMED_REF_SEQ_IDS_TXT}_trimmed_ids.txt"

POSTFIX="${ACCESSIONS//:/_}"

# Raw query sequences from NCBI database
ACC_LIST="$DATA_DIR/accessions_${POSTFIX}.txt"
QUERY_SEQS_GB="$DATA_DIR/query_seqs_${POSTFIX}.gb"
QUERY_SEQS_FASTQ_DIR="$BUILD_DIR/query_seqs_${POSTFIX}"
QUERY_SEQS_MANIFEST_CSV="$BUILD_DIR/query_manifest_${POSTFIX}.csv"
QUERY_SEQS_QZA="$BUILD_DIR/query_seqs_${POSTFIX}.qza"
QUERY_SEQS_TABLE_QZA="$BUILD_DIR/query_table_${POSTFIX}.qza"
QUERY_SEQS_DEREP_QZA="$BUILD_DIR/query_seqs_derep_${POSTFIX}.qza"

# Clustered query sequences
CLUST_PREFIX="gg_13_8_97ref_97clust"
CLUST_DIR="$BUILD_DIR/$CLUST_PREFIX"
CLUST_TABLE_QZA="$CLUST_DIR/${CLUST_PREFIX}_table.qza"
CLUST_TABLE_TSV="$CLUST_DIR/${CLUST_PREFIX}_table.tsv"
CLUST_SEQS_QZA="$CLUST_DIR/${CLUST_PREFIX}_seqs.qza"
CLUST_SEQS_FA="$CLUST_DIR/${CLUST_PREFIX}_seqs.fa"
CLUST_NEW_REFS_QZA="$CLUST_DIR/${CLUST_PREFIX}_new_refs.qza"
CLUST_NEW_REFS_FA="$CLUST_DIR/${CLUST_PREFIX}_new_refs.fa"

# Chimera filtered outputs
FILTERED_PREFIX="${CLUST_PREFIX}_filtered"
FILTERED_DIR="$BUILD_DIR/${FILTERED_PREFIX}"
FILTERED_TABLE_QZA="$FILTERED_DIR/${FILTERED_PREFIX}_table.qza"
FILTERED_TABLE_TSV="$DATA_DIR/${FILTERED_PREFIX}_table.tsv"
FILTERED_SEQS_QZA="$FILTERED_DIR/${FILTERED_PREFIX}_seqs.qza"
FILTERED_SEQS_FA="$DATA_DIR/${FILTERED_PREFIX}_seqs.fa"
FILTERED_NEW_REFS_FA="$DATA_DIR/${FILTERED_PREFIX}_new_refs.fa"
FILTERED_NEW_REFS_ALIGNED_FA="$DATA_DIR/${FILTERED_PREFIX}_new_refs_aligned.fa"

# Fragment insertion output directory
INSERT_PREFIX="${FILTERED_PREFIX}_extended"
INSERT_DIR="$BUILD_DIR/${INSERT_PREFIX}"
REF_TREE_BASE=$(basename "$REF_TREE_PRUNED_NWK")
REF_TREE_BASE="${REF_TREE_BASE%.nwk}"
REF_MODEL_PREFIX="$INSERT_DIR/${REF_TREE_BASE}"
REF_MODEL_STATS="${REF_MODEL_PREFIX}.raxml.bestModel"

# --- Helpers ---

step() {
  echo -e "\n=== Step $1: $2 ==="
}

# --- Main ---

mkdir -p "$BUILD_DIR" "$QUERY_SEQS_FASTQ_DIR" "$CLUST_DIR" "$FILTERED_DIR" "$INSERT_DIR"

# # Prune the reference sequences and tree

# step 1 "Reduce alignment"
# fasta-utils/reduce_alignment.sh \
#   --i-ref-sequences $REF_SEQS_FA \
#   --i-ref-alignment $REF_SEQS_ALIGNED_FA \
#   --o-ref-sequences-reduced $REF_SEQS_REDUCED_FA \
#   --o-ref-alignment-reduced $REF_SEQS_ALIGNED_REDUCED_FA \
#   --threshold $THRESH

# step 2 "Prune reference tree"
# nwk-utils/prune_tree.sh \
#   --i-ref-tree $REF_TREE_NWK \
#   --i-ids-to-prune $TRIMMED_REF_SEQ_IDS_TXT \
#   --o-ref-tree-pruned $REF_TREE_PRUNED_NWK

# # Download GenBank files from NCBI

# step 3 "Build accession list"
# fasta-utils/fetch_accession_list.sh \
#   --i-accession "$ACCESSIONS" \
#   --o-outdir "$DATA_DIR"

# step 4 "Fetch GenBank files"
# fasta-utils/fetch_ncbi_genbank_files.sh \
#   --i-accession-list "$ACC_LIST" \
#   --o-genbank "$QUERY_SEQS_GB" \
#   --batch-size 100

# # Perform open-ref clustering and chimera filtering in QIIME2

# step 5 "Convert GenBank to FASTQ and build manifest"
# python3 qiime2-utils/build_fastq_manifest.py \
#   --i-genbank "$QUERY_SEQS_GB" \
#   --o-outdir "$QUERY_SEQS_FASTQ_DIR" \
#   --o-manifest "$QUERY_SEQS_MANIFEST_CSV"

# step 6 "Import query sequences into QIIME2"
# qiime2-utils/import_seqs_from_manifest.sh \
#   --i-manifest "$QUERY_SEQS_MANIFEST_CSV" \
#   --o-sequences "$QUERY_SEQS_QZA"

# step 7 "Dereplicate sequences"
# qiime2-utils/dereplicate_seqs.sh \
#   --i-sequences "$QUERY_SEQS_QZA" \
#   --o-table "$QUERY_SEQS_TABLE_QZA" \
#   --o-sequences "$QUERY_SEQS_DEREP_QZA"

# step 8 "Open-reference clustering"
# qiime2-utils/open_ref_cluster_seqs.sh \
#   --i-ref-fasta "$REF_SEQS_REDUCED_FA" \
#   --i-table "$QUERY_SEQS_TABLE_QZA" \
#   --i-query-qza "$QUERY_SEQS_DEREP_QZA" \
#   --i-identity 0.97 \
#   --prefix "$CLUST_PREFIX" \
#   --o-outdir "$CLUST_DIR"

# step 9 "Convert table to TSV"
# qiime2-utils/export_table_to_tsv_from_qza.sh \
#   -q "$CLUST_TABLE_QZA" \
#   -t "$CLUST_TABLE_TSV"

# step 10 "Export new reference sequences to FASTA"
# qiime2-utils/export_seqs_to_fasta_from_qza.sh \
#   -q "$CLUST_SEQS_QZA" \
#   -f "$CLUST_SEQS_FA"

# step 11 "Export de-novo clustered sequences to FASTA"
# qiime2-utils/export_seqs_to_fasta_from_qza.sh \
#   -q "$CLUST_NEW_REFS_QZA" \
#   -f "$CLUST_NEW_REFS_FA"

# step 12 "Filter chimeric sequences"
# qiime2-utils/filter_chimeras.sh \
#   --i-sequences "$CLUST_SEQS_QZA" \
#   --i-table "$CLUST_TABLE_QZA" \
#   --i-ref-chimera-fasta "$REF_CHIMERA_FA" \
#   --prefix "$FILTERED_PREFIX" \
#   --o-outdir "$FILTERED_DIR"

# step 13 "Convert chimera-filtered table to TSV"
# qiime2-utils/export_table_to_tsv_from_qza.sh \
#   -q "$FILTERED_TABLE_QZA" \
#   -t "$FILTERED_TABLE_TSV"

# step 14 "Export chimera-filtered new reference sequences to FASTA"
# qiime2-utils/export_seqs_to_fasta_from_qza.sh \
#   -q "$FILTERED_SEQS_QZA" \
#   -f "$FILTERED_SEQS_FA"

# # Pre-process clustered query sequences

# step 15 "Parse out new reference sequences"
# fasta-utils/parse_new_refs.sh \
#   --i-ref-sequences $REF_SEQS_REDUCED_FA \
#   --i-query-sequences $FILTERED_SEQS_FA \
#   --o-new-ref-sequences $FILTERED_NEW_REFS_FA \
#   --o-outdir $DATA_DIR

step 16 "Align new reference sequences"
fasta-utils/align_new_refs.sh

# Compute model statistics with RAxML

step 17 "Compute model statistics"
echo -e "${CYAN}NOTE: This step can take a while on large trees...${NC}"
raxml-utils/compute_model_stats.sh \
  --i-ref-tree $REF_TREE_PRUNED_NWK \
  --i-ref-alignment $REF_SEQS_ALIGNED_REDUCED_FA \
  --o-outdir $INSERT_DIR 

# Perform fragment insertion with SEPP

step 18 "Fragment insertion"
epa-ng-utils/insert_fragments.sh \
  --i-ref-tree $REF_TREE_PRUNED_NWK \
  --i-ref-alignment $REF_SEQS_ALIGNED_REDUCED_FA \
  --i-ref-model-stats $REF_MODEL_STATS \
  --i-query-sequences $FILTERED_NEW_REFS_FA \
  --o-outdir $INSERT_DIR \
  --prefix $INSERT_PREFIX

echo "${CYAN}Chimera-filtered OTU Table:${NC} $FILTERED_TABLE_TSV"
echo "${CYAN}Chimera-filtered OTU Sequences:${NC} $FILTERED_SEQS_FA"
echo "${CYAN}New OTUs Added:${NC} $FILTERED_NEW_REFS_FA"
echo "${CYAN}Extended reference tree in:${NC} $INSERT_DIR"
echo -e "${GREEN}\nPipeline complete.${NC}"