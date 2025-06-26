#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# --- Parse inputs ---

INPUT_QZA=""        # Input QZA filename
OUTPUT_FASTA=""     # Output FASTA filename

usage() {
  echo -e "${CYAN}Usage:${NC} $0 [-q input.qza | --qza input.qza] [-f output.fasta | --fasta output.fasta]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -q|--qza) 
      INPUT_QZA="$2" 
      shift 2 
      ;;
    -f|--fasta|--fa) 
      OUTPUT_FASTA="$2" 
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

if [[ -z "$INPUT_QZA" ]] || [[ -z "$OUTPUT_FASTA" ]]; then
  echo -e "${RED}Error:${NC} Both -q and -f are required." >&2
  usage
fi

# --- Perform checks --

if [[ ! -f "$INPUT_QZA" ]]; then
  echo -e "${RED}Error:${NC} Input file not found: $INPUT_QZA" >&2
  exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
  echo -e "${RED}Error:${NC} 'qiime' command not found. Please install QIIME 2 and ensure it is in your PATH." >&2
  exit 1
fi

# --- Main ---

EXPORT_DIR=$(mktemp -d)

echo -e "${CYAN}Exporting QZA artifact...${NC}"
qiime tools export \
  --input-path "$INPUT_QZA" \
  --output-path "$EXPORT_DIR"

FASTA_FILE=$(find "$EXPORT_DIR" -name '*.fasta' | head -n 1)

if [[ -f "$FASTA_FILE" ]]; then
  cp "$FASTA_FILE" "$OUTPUT_FASTA"
else
  echo -e "${RED}Error:${NC} No FASTA file found inside exported QZA artifact." >&2
  rm -rf "$EXPORT_DIR"
  exit 1
fi

echo -e "${CYAN}FASTA file saved to: $OUTPUT_FASTA${NC}"
echo -e "${CYAN}Cleaning up temporary directory...${NC}"
rm -rf "$EXPORT_DIR"
echo -e "${GREEN}Done converting QZA to FASTA.${NC}"
