#!/usr/bin/env bash

set -euo pipefail

buildCentroidComparisonReport.pl --oldBestReps $oldBestReps \\
                                  --newBestReps $newBestReps \\
                                  --oldStats $oldStats \\
                                  --newStats $newStats \\
                                  --centroidsTsv $centroidsTsv \\
                                  --groupsFile $groupsFile \\
                                  --outputFile centroidComparison.tsv
