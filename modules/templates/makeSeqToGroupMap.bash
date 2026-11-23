#!/usr/bin/env bash

set -euo pipefail

makeSeqToGroupMap.pl --groupFile $groupFile --outputFile seqToGroup.tsv
