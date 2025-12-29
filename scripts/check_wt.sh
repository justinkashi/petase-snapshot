# usage: ./check_wt.sh "A214H/I168R/..." "SEQUENCESTRING"
# checks that the WT (first letter) matches the input sequence at each position
# goes through ALL mutations; exits nonzero if any fail

set -euo pipefail

MUTS="$1"
SEQ="$2"
len=${#SEQ}

fails=0

IFS='/' read -ra arr <<< "$MUTS"
for mut in "${arr[@]}"; do
  if [[ ! "$mut" =~ ^([A-Z])([0-9]+)([A-Z])$ ]]; then
    echo "Invalid mutation format: $mut" >&2
    fails=$((fails+1))
    continue
  fi

  wt="${BASH_REMATCH[1]}"
  pos="${BASH_REMATCH[2]}"
  aa="${BASH_REMATCH[3]}"

  if (( pos < 1 || pos > len )); then
    echo "Position $pos out of range (seq length $len) for $mut" >&2
    fails=$((fails+1))
    continue
  fi

  seq_aa="${SEQ:pos-1:1}"

  if [[ "$seq_aa" == "$wt" ]]; then
    echo "OK  $mut → sequence[$pos]=$seq_aa (WT matches)"
  else
    echo "FAIL $mut → sequence[$pos]=$seq_aa (expected WT=$wt)" >&2
    fails=$((fails+1))
  fi
done

(( fails == 0 )) || exit 1