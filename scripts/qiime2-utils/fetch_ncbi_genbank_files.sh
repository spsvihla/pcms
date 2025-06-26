#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

ACCESSION_LIST=""                # Accession list (required)
GENBANK="query_seqs.gb"          # Output GenBank filename 
BATCH_SIZE=100                   # Batch size when fetching from NCBI

usage() {
  echo "${CYAN}Usage:${NC} $0 --i-accession-list accession_list.txt [--o-genbank output_genbank_file] [-b batch_size | --batch-size batch_size]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-accession-list)
      ACCESSION_LIST="$2"
      shift 2
      ;;
    --o-genbank)
      GENBANK="$2"
      shift 2
      ;;
    -b|--batch-size)
      BATCH_SIZE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo -e "${RED}Invalid option: $1${NC}" >&2
      usage
      ;;
  esac
done

if [[ -z "$ACCESSION_LIST" ]]; then
  echo -e "${RED}Error:${NC} --i-accession-list argument is required." >&2
  usage
fi

# --- Perform checks ---

for cmd in epost efetch split realpath; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo -e "${RED}Error:${NC} Required command '$cmd' not found in PATH." >&2
    exit 1
  fi
done

# --- Main ---

ACCESSION_LIST="$(realpath "$ACCESSION_LIST")"
GENBANK="$(realpath -m "$GENBANK")" # -m allows non-existent target

mkdir -p "$(dirname "$GENBANK")"
if [[ -f "$GENBANK" ]]; then
  echo -e "${CYAN}Output file '$GENBANK' already exists. Skipping fetch.${NC}"
  exit 0
fi

TMP_DIR="$(dirname "$GENBANK")/batch_tmp"
if [[ -d "$TMP_DIR" ]]; then
  echo -e "${CYAN}Cleaning up existing temporary directory:${NC} $TMP_DIR"
  rm -rf "$TMP_DIR"
fi
mkdir "$TMP_DIR"

echo -e "${CYAN}Splitting accession list into batches of $BATCH_SIZE in $TMP_DIR$...${NC}"
split -l "$BATCH_SIZE" "$ACCESSION_LIST" "$TMP_DIR/batch_"

num_batches=$(ls "$TMP_DIR"/batch_* | wc -l)
echo -e "${CYAN}Starting fetch loop over $num_batches batch(es)...${NC}"

count=0
for batch_file in "$TMP_DIR"/batch_*; do
  count=$((count + 1))
  echo -e "${CYAN}[$count/$num_batches] Fetching batch: $(basename "$batch_file")${NC}"

  if ! cat "$batch_file" | epost -db nucleotide | efetch -format gbwithparts >> "$GENBANK"; then
    echo -e "${RED}Warning:${NC} Failed to fetch batch: $(basename "$batch_file")" >&2
  fi

  sleep 1
done

echo -e "${CYAN}Fetched sequences saved to:${NC} $GENBANK"
echo -e "${CYAN}Cleaning up temporary directory...${NC}"
rm -rf "$TMP_DIR"
echo -e "${GREEN}Done fetching GenBank sequence files.${NC}"
