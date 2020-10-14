#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# setup a list of the different test_sets to run
test_sets=(100_species 10_species 500_species 5638_species 1000_species)

for test_set in "${test_sets[@]}"; do

  # benchmark malt building db performance test for dbs of various size
  /usr/bin/time -v malt-build \
    -t "$MAX_CPU" \
    -s DNA \
    -i db_mutlifasta_inputs/db_input_"${test_set}".fasta \
    -d index_new_"${test_set}"/ \
    -a2taxonomy mapping_files/megan-map-Jul2020-2.db

done
