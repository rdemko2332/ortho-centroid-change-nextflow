#!/usr/bin/env bash

set -euo pipefail

for f in OG*.fasta; do
  [ -f "\$f" ] || continue
  basename "\$f" .fasta
done > batchGroupIds.txt

cat OG*.fasta > batchMembers.fasta

samtools faidx $newBestRepsFull
samtools faidx $newBestRepsFull -r batchGroupIds.txt > batchCentroids.fasta

diamond makedb --in batchCentroids.fasta --db batchCentroids

diamond blastp -d batchCentroids.dmnd -q batchMembers.fasta -o batchHits.tsv \\
  -e $distanceEvalue -k 0 -f 6 $diamondOutputFields \\
  --comp-based-stats 0 --very-sensitive --ignore-warnings
