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
QUERY_SEQS_FA=""
OUTDIR="."
PREFIX=""

usage() {
  echo "Usage: $0 --i-ref-tree tree.nwk --i-ref-alignment ref_seqs_aligned.fasta --i-ref-model-stats ref_model_stats.raxml.bestModel --i-query-sequences query_seqs.fasta [--i-outdir outdir/] [-p prefix | --prefix prefix]"
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
    --i-ref-model-stats)
      [[ $# -lt 2 ]] && usage
      REF_MODEL_STATS="$2"
      shift 2 
      ;;
    --i-query-sequences) 
      [[ $# -lt 2 ]] && usage
      QUERY_SEQS_FA="$2" 
      shift 2
      ;;
    --o-outdir)
      [[ $# -lt 2 ]] && usage
      OUTDIR="$2"
      shift 2
      ;;
    -p|--prefix)
      [[ $# -lt 2 ]] && usage
      PREFIX="$2"
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

if [[ -z "$REF_TREE_NWK" || -z "$REF_SEQS_ALIGNED_FA" || -z "$REF_MODEL_STATS" || -z "$QUERY_SEQS_FA" ]]; then
  echo -e "${RED}Error:${NC} --i-ref-tree tree.nwk --i-ref-alignment ref_seqs_aligned.fasta --i-ref-model-stats ref_model_stats.raxml.bestModel and --i-query-sequences query_seqs.fasta are required." >&2
  usage
fi

# --- Perform checks ---

command -v epa-ng >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'epa-ng' command not found in PATH. Please install EPA-NG and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

# NUM_CPUS=$(nproc)
NUM_CPUS=12

echo -e "${CYAN}Constructing updated tree...${NC}"
epa-ng \
  --tree $REF_TREE_NWK \
  --ref-msa $REF_SEQS_ALIGNED_FA \
  --query $QUERY_SEQS_FA \
  --model $REF_MODEL_STATS \
  --outdir $OUTDIR \
  --preserve-rooting on \
  --threads $NUM_CPUS \
  --redo

echo -e "${CYAN}Extended tree located in:${NC} $OUTDIR"
echo -e "${GREEN}Done inserting fragments.${NC}"
