#!/usr/bin/env bash

set -euo pipefail

applyReassignments.pl --groupFile $groupFile \\
                       --reassignments $reassignments \\
                       --outputFile newGroupFile.txt \\
                       --dissolvedSingletons dissolvedSingletons.tsv
