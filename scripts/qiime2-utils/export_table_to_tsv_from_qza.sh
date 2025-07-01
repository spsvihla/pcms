#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# --- Parse inputs ---

INPUT_QZA=""      # Input QZA filename
OUTPUT_TSV=""     # Output TSV filename

usage() {
  echo -e "${CYAN}Usage:${NC} $0 [-q input.qza | --qza input.qza] [-t output.tsv | --tsv output.tsv]"
  exit 1
}

# Parse options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -q|--qza) 
      [[ $# -lt 2 ]] && usage
      INPUT_QZA="$2" 
      shift 2
      ;;
    -t|--tsv) 
      [[ $# -lt 2 ]] && usage
      OUTPUT_TSV="$2"
      shift 2
      ;;
    -h|--help) 
      usage 
      ;;
    *) 
      echo -e "${RED}Invalid option: -$1${NC}" >&2; 
      usage 
      ;;
  esac
done

if [[ -z "$INPUT_QZA" ]] || [[ -z "$OUTPUT_TSV" ]]; then
  echo -e "${RED}Error:${NC} Both -q and -t are required." >&2
  usage
fi

# --- Perform checks ---

if [[ ! -f "$INPUT_QZA" ]]; then
  echo -e "${RED}Error:${NC} Input file not found: $INPUT_QZA" >&2
  exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
  echo -e "${RED}Error:${NC} 'qiime' command not found in PATH. Please install QIIME 2 and ensure it is in your PATH." >&2
  exit 1
fi

if ! command -v biom >/dev/null 2>&1; then
  echo -e "${RED}Error:${NC} 'biom' command not found in PATH. Install with 'pip install biom-format'." >&2
  exit 1
fi

# --- Main ---

EXPORT_DIR=$(mktemp -d)

echo -e "${CYAN}Exporting feature table artifact...${NC}"
qiime tools export \
  --input-path "$INPUT_QZA" \
  --output-path "$EXPORT_DIR"

BIOM_FILE="$EXPORT_DIR/feature-table.biom"

if [[ ! -f "$BIOM_FILE" ]]; then
  echo -e "${RED}Error:${NC} Exported BIOM file not found: $BIOM_FILE" >&2
  exit 1
fi

echo -e "${CYAN}Converting BIOM to TSV format...${NC}"
biom convert \
  -i "$BIOM_FILE" \
  -o "$OUTPUT_TSV" \
  --to-tsv

echo -e "${CYAN}TSV file saved to:${NC} $OUTPUT_TSV"
echo -e "${CYAN}Cleaning up temporary directory...${NC}"
rm -rf "$EXPORT_DIR"
echo -e "${GREEN}Done converting QZA to TSV${NC}"
