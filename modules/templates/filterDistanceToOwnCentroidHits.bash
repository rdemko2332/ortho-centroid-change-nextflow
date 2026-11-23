#!/usr/bin/env bash

set -euo pipefail

filterDistanceToOwnCentroidHits.pl --diamondOut $batchHits \\
                                    --seqToGroup $seqToGroup \\
                                    --bestReps $newBestRepsFullTxt \\
                                    --outputFile batchDistances.tsv
