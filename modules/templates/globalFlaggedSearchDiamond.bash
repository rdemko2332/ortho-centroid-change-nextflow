#!/usr/bin/env bash

set -euo pipefail

diamond makedb --in $newBestRepsFull --db newdb

diamond blastp -d newdb.dmnd -q $flaggedSubset -o flaggedHits.tsv \\
  -e $globalSearchEvalue -k 1 -f 6 $diamondOutputFields \\
  --comp-based-stats 0 --no-self-hits --very-sensitive --ignore-warnings
