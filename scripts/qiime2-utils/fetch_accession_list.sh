#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

ACCESSION=""         # Accession range or list input string (required)
OUTPUT_DIR="."       # Output directory default is current directory

usage() {
  echo -e "${CYAN}Usage:${NC} $0 --i-accession accession_range_or_list [--o-directory output_dir]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-accession)
      ACCESSION="$2"
      shift 2
      ;;
    --o-outdir)
      OUTPUT_DIR="$2"
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

if [[ -z "$ACCESSION" ]]; then
  echo -e "${RED}Error:${NC} --i-accession argument is required." >&2
  usage
fi

# --- Perform checks ---

for cmd in esearch efetch; do
  if ! command -v "$cmd" &>/dev/null; then
    echo -e "${RED}Error:${NC} '$cmd' not found in PATH. Please install NCBI EDirect." >&2
    exit 1
  fi
done

# --- Main ---

mkdir -p "$OUTPUT_DIR"

ACCESSION_LIST_BASENAME=$(echo "$ACCESSION" | sed 's/[:,]/_/g')
ACCESSION_LIST="$OUTPUT_DIR/accessions_${ACCESSION_LIST_BASENAME}.txt"

if [[ -f "$ACCESSION_LIST" ]]; then
  echo -e "${CYAN}Output file '$ACCESSION_LIST' already exists. Skipping fetch.${NC}"
  exit 0
fi

echo -e "${CYAN}Downloading accession list for:${NC} \"$ACCESSION\""

esearch -db nucleotide -query "${ACCESSION}[Accession]" | efetch -format acc > "$ACCESSION_LIST"

NUM_ACC=$(wc -l < "$ACCESSION_LIST" | tr -d '[:space:]')
echo -e "${CYAN}Downloaded ${NUM_ACC} accession(s) to${NC} ${ACCESSION_LIST}"
echo -e "${GREEN}Done fetching accession list.${NC}"
