#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# setup a list of the different test_sets to run
test_sets=(100_species 10_species 1000_species 500_species 5638_species)

for test_set in "${test_sets[@]}"; do
  # benchmark kraken building db performance test for 100 species
  /usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_"${test_set}"

  /usr/bin/time -v kraken2-build \
  --add-to-library db_mutlifasta_inputs/db_input_"${test_set}".fasta \
  --db db_kraken_"${test_set}" \
  --threads "$MAX_CPU"

  /usr/bin/time -v kraken2-build --build --db db_kraken_"${test_set}" --threads "$MAX_CPU"

  # run bracken building db performance test for 100 species
  /usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_"${test_set}"/

done