#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

test_sets=(input_10K input_100K input_100M input_1M input_10M)

for test_set in "${test_sets[@]}"; do

  mkdir -p sigma_outputs/"${test_set}"

  # benchmark kraken sample analysis performance test for samples of various sizes
  /usr/bin/time -v kraken2 \
    --db db_kraken_5638_species \
    --gzip-compressed \
    --use-names \
    --report "${test_set}".report \
    --unclassified-out "${test_set}"_unclassified.out \
    --classified-out "${test_set}"_classified.out \
    --threads "$MAX_CPU" \
    --output "${test_set}".out ./inputs/"${test_set}"_reads.fastq.gz; \
    Bracken-2.5/bracken -d db_kraken_5638_species \
    -i "${test_set}".report \
    -o "${test_set}".bracken \
    --threads "$MAX_CPU"

done