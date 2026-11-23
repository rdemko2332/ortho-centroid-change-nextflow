#!/usr/bin/env bash

set -euo pipefail

# Groups present in groupsFile but with zero rows in blastFile at all.
grep -oP '^\\S+(?=:)' $groupsFile | sort -u > allGroups.txt
cut -f1 $blastFile | sort -u > groupsWithHits.txt
comm -23 allGroups.txt groupsWithHits.txt > missingGroups.txt

calculateCentroidDistanceStats.pl --blastFile $blastFile \\
                                   --groupsFile $groupsFile \\
                                   --missingGroups missingGroups.txt \\
                                   --outputFile distance_stats.txt
