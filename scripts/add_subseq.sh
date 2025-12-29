# usage: ./subseq_add.sh "SUBSEQUENCE" "SEQUENCE"
# checks whether SUBSEQUENCE occurs in SEQUENCE (exact match)

set -euo pipefail

SUB="$1"
SEQ="$2"

if [[ -z "$SUB" || -z "$SEQ" ]]; then
  echo "Empty subsequence or sequence" >&2
  exit 1
fi

if [[ "$SEQ" == *"$SUB"* ]]; then
  pos=$(expr index "$SEQ" "$SUB")
  echo "FOUND at position $pos"
else
  echo "NOT FOUND"
  exit 1
fi
