#!/usr/bin/env bash

set -euo pipefail

mkdir -p touchedGroupFastas

splitTouchedGroupFastas.pl --groupFile $groupFile \\
                           --touchedGroups $touchedGroups \\
                           --proteome $fullProteome \\
                           --outputDir touchedGroupFastas
