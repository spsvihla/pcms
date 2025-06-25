#!/usr/bin/env bash

cd ./open_ref_cluster_from_genbank/

export DATASET_DIR=$DATA/guerrero_negro
export BUILD_DIR=$DATASET_DIR/build
mkdir -p "$BUILD_DIR"

# Step 1: Build accession list for the desired range
./build_accession_list.sh \
  -a "JN427016:JN539989" \
  -o "$DATASET_DIR"

# Step 2: Download GenBank files for accessions
./fetch_ncbi_genbank_files.sh \
  -i "$DATASET_DIR/accessions_JN427016_JN539989.txt" \
  -o "$DATASET_DIR/query_seqs_JN427016_JN539989.gb"

export REF_SEQS=$DATA/greengenes/gg_13_8_otus/rep_set/97_otus.fasta
export QUERY_FASTQ_DIR=$BUILD_DIR/query_seqs_JN427016_JN539989_fastq/

# Step 3: Convert GenBank file to FASTA using Python script
python3 build_seq_manifest.py \
  -i "$DATASET_DIR/query_seqs_JN427016_JN539989.gb" \
  -o "$QUERY_FASTQ_DIR" \
  -m "$BUILD_DIR/query_manifest_JN427016_JN539989.csv"

# Step 4: Import query sequences into QIIME2 artifact
./import_query_seqs_from_manifest.sh \
  -m "$BUILD_DIR/query_manifest_JN427016_JN539989.csv" \
  -o "$BUILD_DIR/query_seqs_JN427016_JN539989.qza" 

# Step 5: Dereplicate sequences to generate feature table and representative sequences
./dereplicate_seqs.sh \
  -i "$BUILD_DIR/query_seqs_JN427016_JN539989.qza" \
  -t "$BUILD_DIR/query_table_JN427016_JN539989.qza" \
  -s "$BUILD_DIR/query_seqs_derep_JN427016_JN539989.qza"

# Step 6: Perform open-reference clustering against Greengenes reference
./open_ref_cluster.sh \
  -r "$REF_SEQS" \
  -t "$BUILD_DIR/query_table_JN427016_JN539989.qza" \
  -q "$BUILD_DIR/query_seqs_derep_JN427016_JN539989.qza" \
  -i 0.97 \
  -p "gg_13_8_97ref_97clust_JN427016_JN539989" \
  -o "$BUILD_DIR"

# Step 7: Convert .qza to .tsv
./convert_qza_table2csv.sh \
  -i "$BUILD_DIR/gg_13_8_97ref_97clust_JN427016_JN539989_table.qza" \
  -o "$DATASET_DIR/gg_13_8_97ref_97clust_JN427016_JN539989_table.tsv"