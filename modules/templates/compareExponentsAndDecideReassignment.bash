#!/usr/bin/env bash

set -euo pipefail

compareExponentsAndDecideReassignment.pl --flaggedProteins $flaggedProteins \\
                                          --flaggedHits $flaggedHits \\
                                          --exponentImprovement $exponentImprovement \\
                                          --outputFile reassignments.tsv
