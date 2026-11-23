#!/usr/bin/env bash

set -euo pipefail

mkdir -p exact large

for f in OG*.fasta; do
  [ -f "\$f" ] || continue
  SEQ_COUNT=\$(grep -c ">" "\$f")

  if [ "\$SEQ_COUNT" -le 2 ]; then
    # Untouched -- step 1 only recomputes groups with more than 2 members.
    continue
  elif [ "\$SEQ_COUNT" -le $exactCentroidLimit ]; then
    cp "\$f" exact/
  else
    cp "\$f" large/
  fi
done
