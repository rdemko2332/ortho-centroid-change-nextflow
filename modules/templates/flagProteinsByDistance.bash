#!/usr/bin/env bash

set -euo pipefail

flagProteinsByDistance.pl --bestRepFile $newBestRepsFullTxt \\
                           --groupsFile $groupsFile \\
                           --blastFile $blastFile \\
                           --threshold $threshold \\
                           --outputFile flaggedProteins.tsv
