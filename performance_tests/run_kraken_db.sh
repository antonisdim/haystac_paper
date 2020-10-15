#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
test_set=${args[0]}

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# benchmark kraken building db performance test for 100 species
kraken2-build --download-taxonomy --db db_kraken_"${test_set}"

kraken2-build \
  --add-to-library db_mutlifasta_inputs/db_input_"${test_set}".fasta \
  --db db_kraken_"${test_set}" \
  --threads "$MAX_CPU"

kraken2-build --build --db db_kraken_"${test_set}" --threads "$MAX_CPU"

# run bracken building db performance test for 100 species
Bracken-2.5/bracken-build -d db_kraken_"${test_set}"/

