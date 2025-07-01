#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

REF_TREE_NWK=""
REF_SEQS_ALIGNED_FA=""
OUTDIR="."

usage() {
  echo "Usage: $0 --i-ref-tree tree.nwk --i-ref-alignment ref_seqs_aligned.fasta [--o-outdir outdir]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-ref-tree) 
      [[ $# -lt 2 ]] && usage
      REF_TREE_NWK="$2" 
      shift 2
      ;;
    --i-ref-alignment) 
      [[ $# -lt 2 ]] && usage
      REF_SEQS_ALIGNED_FA="$2" 
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
    *)
      echo -e "${RED}Invalid option: -$1${NC}" >&2
      usage
      ;;
  esac
done

if [[ -z "$REF_TREE_NWK" || -z "$REF_SEQS_ALIGNED_FA" ]]; then
  echo -e "${RED}Error:${NC} --i-ref-tree tree.nwk and --i-ref-alignment ref_seqs_aligned.fasta are required." >&2
  usage
fi

# --- Perform checks ---

command -v raxml-ng >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'raxml-ng' command not found in PATH. Please install RAxML-NG and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

# NUM_CPUS=$(nproc)
NUM_CPUS=12

REF_TREE_BASE=$(basename "$REF_TREE_NWK")
REF_TREE_BASE="${REF_TREE_BASE%.nwk}"
REF_MODEL_PREFIX="$OUTDIR/${REF_TREE_BASE}"
REF_MODEL_STATS="${REF_MODEL_PREFIX}.raxml.bestModel"

WORKING_DIR=$(pwd)
mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ ! -f "$REF_MODEL_STATS" ]]; then
  echo -e "${CYAN}Computing tree statistics with RAxML...${NC}"
  raxml-ng --evaluate \
  --msa "$REF_SEQS_ALIGNED_FA" \
  --msa-format FASTA \
  --model GTR+G \
  --tree "$REF_TREE_NWK" \
  --data-type DNA \
  --threads "$NUM_CPUS" \
  --prefix "$REF_MODEL_PREFIX" \
  --pat-comp on \
  --tip-inner on \
  --simd avx2 \
  --opt-model on \
  --opt-branches on \
  --redo
else
  echo -e "${CYAN}RAxML stats file already exists at $OUTDIR/${REF_MODEL_STATS}, skipping computation.${NC}"
fi

cd $WORKING_DIR

echo -e "${GREEN}Done computing model statistics.${NC}"