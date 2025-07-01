#!/usr/bin/env bash
set -euo pipefail

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Parse inputs ---

REF_TREE_NWK=""
TRIMMED_REF_SEQ_IDS_TXT=""
REF_TREE_PRUNED_NWK=""

usage() {
  echo "Usage: $0 --i-ref-tree ref_tree.nwk --i-ids-to-prune ids.txt [--o-ref-tree-pruned ref_tree_pruned.nwk]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i-ref-tree)
      [[ $# -lt 2 ]] && usage
      REF_TREE_NWK="$2" 
      shift 2
      ;;
    --i-ids-to-prune) 
      [[ $# -lt 2 ]] && usage
      TRIMMED_REF_SEQ_IDS_TXT="$2" 
      shift 2
      ;;
    --o-ref-tree-pruned)
      [[ $# -lt 2 ]] && usage
      REF_TREE_PRUNED_NWK="$2" 
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

if [[ -z "$REF_TREE_NWK" || -z "$TRIMMED_REF_SEQ_IDS_TXT" ]]; then
  echo -e "${RED}Error:${NC} --i-ref-tree and --i-ids-to-prune are required." >&2
  usage
fi

if [[ -z "$REF_TREE_PRUNED_NWK" ]]; then
    REF_TREE_PRUNED_NWK="${REF_TREE_NWK%.nwk}"
    REF_TREE_PRUNED_NWK="${REF_TREE_PRUNED_NWK}_pruned.nwk"
fi

# --- Perform checks ---

command -v nw_prune >/dev/null 2>&1 || {
  echo -e "${RED}Error:${NC} 'nw_prune' command not found in PATH. Please install 'nw_prune' and ensure it is accessible." >&2
  exit 1
}

# --- Main ---

echo -e "${CYAN}Pruning reference tree.${NC}"

if [[ ! -s "$TRIMMED_REF_SEQ_IDS_TXT" ]]; then
  echo -e "${CYAN}No sequence IDs to prune. Copying original tree to ${REF_TREE_PRUNED_NWK}${NC}"
  cp "$REF_TREE_NWK" "$REF_TREE_PRUNED_NWK"
else
  echo -e "${CYAN}Pruning tree ${REF_TREE_NWK} by removing sequences listed in ${TRIMMED_REF_SEQ_IDS_TXT}${NC}"
  nw_prune -f "$REF_TREE_NWK" "$TRIMMED_REF_SEQ_IDS_TXT" > "$REF_TREE_PRUNED_NWK"
fi

echo -e "${CYAN}Pruned tree saved to:${NC} ${REF_TREE_PRUNED_NWK}"
echo -e "${GREEN}Done trimming alignment sequences.${NC}"