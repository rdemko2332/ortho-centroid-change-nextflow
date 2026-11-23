#!/usr/bin/env bash

set -euo pipefail

cut -f1 $flaggedProteins | sort -u > flaggedIds.txt

samtools faidx $fullProteome
samtools faidx $fullProteome -r flaggedIds.txt > flagged.fasta
