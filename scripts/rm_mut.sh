# usage: ./remove_subseq.sh "SUBSEQUENCE" "SEQUENCE"
# removes the FIRST occurrence of SUBSEQUENCE from SEQUENCE
# intended for signal peptide removal

set -euo pipefail

SUB="$1"
SEQ="$2"

if [[ -z "$SUB" || -z "$SEQ" ]]; then
  echo "Empty subsequence or sequence" >&2
  exit 1
fi

if [[ "$SEQ" != *"$SUB"* ]]; then
  echo "Subsequence not found" >&2
  exit 1
fi

# split at first occurrence
prefix="${SEQ%%"$SUB"*}"
suffix="${SEQ#*"$SUB"}"

# strict signal-peptide check: must be N-terminal
if [[ -n "$prefix" ]]; then
  echo "Error: subsequence not at N-terminus" >&2
  exit 1
fi

echo "$suffix"
