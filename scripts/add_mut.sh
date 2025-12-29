# usage: ./apply_mut.sh "F43D/S142H" "SEQUENCESTRING"
# applies mutations and outputs the mutated sequence
# enforces that WT residue matches the input sequence

set -euo pipefail

MUTS="$1"
SEQ="$2"

len=${#SEQ}
mut_seq="$SEQ"

IFS='/' read -ra arr <<< "$MUTS"
for mut in "${arr[@]}"; do
  if [[ ! "$mut" =~ ^([A-Z])([0-9]+)([A-Z])$ ]]; then
    echo "Invalid mutation format: $mut" >&2
    exit 1
  fi

  wt="${BASH_REMATCH[1]}"
  pos="${BASH_REMATCH[2]}"
  aa="${BASH_REMATCH[3]}"

  if (( pos < 1 || pos > len )); then
    echo "Position $pos out of range (seq length $len)" >&2
    exit 1
  fi

  idx=$((pos - 1))
  current="${mut_seq:idx:1}"

  if [[ "$current" != "$wt" ]]; then
    echo "WT mismatch at $pos: expected $wt, found $current" >&2
    exit 1
  fi

  mut_seq="${mut_seq:0:idx}$aa${mut_seq:idx+1}"
done

echo "$mut_seq"