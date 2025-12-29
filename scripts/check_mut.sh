#!/usr/bin/env bash
# usage: ./check_mut.sh "F43D/S142H" "SEQUENCESTRING"
# checks that the *mutant* residue matches the sequence at each position

set -euo pipefail

MUTS="$1"
SEQ="$2"

len=${#SEQ}

IFS='/' read -ra arr <<< "$MUTS"
for mut in "${arr[@]}"; do
  if [[ ! "$mut" =~ ^([A-Z])([0-9]+)([A-Z])$ ]]; then
    echo "Invalid mutation format: $mut"
    exit 1
  fi

  wt="${BASH_REMATCH[1]}"
  pos="${BASH_REMATCH[2]}"
  aa="${BASH_REMATCH[3]}"

  if (( pos < 1 || pos > len )); then
    echo "Position $pos out of range (seq length $len)"
    exit 1
  fi

  seq_aa="${SEQ:pos-1:1}"

  if [[ "$seq_aa" == "$aa" ]]; then
    echo "OK  $mut → sequence[$pos]=$seq_aa"
  else
    echo "FAIL $mut → sequence[$pos]=$seq_aa (expected $aa)"
  fi
done
