#!/usr/bin/env bash

set -euo pipefail

touch centroids_batch.tsv
touch newCentroids.fasta

for f in OG*.fasta; do
  [ -f "\$f" ] || continue
  groupId=\$(basename "\$f" .fasta)

  python3 $approxCentroidScript "\$f" \\
    --group-id "\$groupId" \\
    --exact-limit $exactCentroidLimit \\
    --threads $task.cpus \\
    --diamond-args '$diamondArgs' \\
    --emit-centroid-seq newCentroids.fasta \\
    | tail -n +2 >> centroids_batch.tsv
done
