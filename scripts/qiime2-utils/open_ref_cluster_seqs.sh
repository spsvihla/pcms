#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

IDENTITY=0.97
PREFIX="cluster"
OUTDIR="."

usage() {
  echo "Usage: $0 "
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-ref-fasta) 
      [[ $# -lt 2 ]] && usage
      REF_FASTA="$2" 
      shift 2
      ;;
    --i-table) 
      [[ $# -lt 2 ]] && usage
      QUERY_TABLE="$2" 
      shift 2
      ;;
    --i-query-qza) 
      [[ $# -lt 2 ]] && usage
      QUERY_QZA="$2" 
      shift 2
      ;;
    --i-identity) 
      [[ $# -lt 2 ]] && usage
      IDENTITY="$2" 
      shift 2
      ;;
    -p|--prefix) 
      [[ $# -lt 2 ]] && usage
      PREFIX="$2" 
      shift 2
      ;;
    --o-outdir) 
      [[ $# -lt 2 ]] && usage
      OUTDIR="$2" 
      shift 2
      ;;
    -h|--help) 
      usage 
      ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

if [[ -z "${REF_FASTA:-}" ]] || [[ -z "${QUERY_QZA:-}" ]] || [[ -z "${QUERY_TABLE:-}" ]]; then
  echo -e "${RED}Error:${NC} --ref-fasta --i-table and --i-query-qza are required." >&2
  usage
fi

# --- Perform checks ---

command -v qiime >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'qiime' command not found in PATH. Please install QIIME 2 and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

mkdir -p "${OUTDIR}"

# Derive reference QZA file name from FASTA base name
REF_BASE=$(basename "$REF_FASTA")
REF_BASE="${REF_BASE%.fasta}"
REF_BASE="${REF_BASE%.fa}"
REF_QZA="${OUTDIR}/${REF_BASE}.qza"

# Define output artifact paths
TABLE_QZA="${OUTDIR}/${PREFIX}_table.qza"
SEQS_QZA="${OUTDIR}/${PREFIX}_seqs.qza"
NEW_REFS_QZA="${OUTDIR}/${PREFIX}_new_refs.qza"
# TABLE_QZV="${OUTDIR}/${PREFIX}_table.qzv"

echo -e "${CYAN}Importing reference sequences from '$REF_FASTA' into QIIME artifact $REF_QZA...${NC}"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$REF_FASTA" \
  --output-path "$REF_QZA" \
  --input-format DNAFASTAFormat

echo -e "${CYAN}Running open-reference clustering (identity threshold: $IDENTITY)...${NC}"
qiime vsearch cluster-features-open-reference \
  --i-table "$QUERY_TABLE" \
  --i-sequences "$QUERY_QZA" \
  --i-reference-sequences "$REF_QZA" \
  --p-perc-identity "$IDENTITY" \
  --o-clustered-table "$TABLE_QZA" \
  --o-clustered-sequences "$SEQS_QZA" \
  --o-new-reference-sequences "$NEW_REFS_QZA"

# echo -e "${CYAN}Generating feature table summary visualization${NC} $TABLE_QZV${CYAN}...${NC}"
# qiime feature-table summarize \
#   --i-table "$TABLE_QZA" \
#   --o-visualization "$TABLE_QZV"

echo -e "${CYAN}Results saved in: ${OUTDIR}${NC}"
echo -e "${GREEN}OTU picking complete.${NC}"
