#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

args=("$@")
test_set=${args[0]}

mkdir -p index_new_"${test_set}"

# benchmark malt building db performance test for dbs of various size
malt-build \
  -t "$MAX_CPU" \
  -s DNA \
  -i db_mutlifasta_inputs/db_input_"${test_set}".fasta \
  -d index_new_"${test_set}"/ \
  -a2taxonomy mapping_files/megan-map-Jul2020-2.db
