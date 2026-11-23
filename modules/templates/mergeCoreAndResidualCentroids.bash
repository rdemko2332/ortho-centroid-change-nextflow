#!/usr/bin/env bash

set -euo pipefail

# Combined group -> rep map: core/peripheral (already merged with untouched
# groups in step 1) + residual (unchanged, per spec step 1 doesn't touch it).
cat $newCoreBestReps $residualBestReps > newBestRepsFull.txt

cut -f2 newBestRepsFull.txt > repSeqIds.txt

samtools faidx $fullProteome
samtools faidx $fullProteome -r repSeqIds.txt \\
  | makeBestRepresentativesFasta.pl --bestReps newBestRepsFull.txt --outputFile newBestRepsFull.fasta
