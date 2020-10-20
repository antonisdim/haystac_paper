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

mkdir -p sigma_outputs/"${test_set}"

# benchmark kraken sample analysis performance test for samples of various sizes
kraken2 \
  --db db_kraken_5638_species \
  --gzip-compressed \
  --use-names \
  --report "${test_set}".report \
  --unclassified-out "${test_set}"_unclassified.out \
  --classified-out "${test_set}"_classified.out \
  --threads "$MAX_CPU" \
  --output "${test_set}".out ./inputs/"${test_set}"_reads.fastq.gz

Bracken-2.5/bracken -d db_kraken_5638_species \
  -i "${test_set}".report \
  -o "${test_set}".bracken \
  --threads "$MAX_CPU"
