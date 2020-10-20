#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
sample_input=${args[0]}
db_input=${args[1]}

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# benchmark kraken sample analysis performance test for samples of various sizes
kraken2 \
  --db db_kraken_"${db_input}"_species \
  --gzip-compressed \
  --use-names \
  --report "${sample_input}"_"${db_input}"_species.report \
  --unclassified-out "${sample_input}"_"${db_input}"_species_unclassified.out \
  --classified-out "${sample_input}"_"${db_input}"_species_classified.out \
  --threads "$MAX_CPU" \
  --output "${sample_input}"_"${db_input}"_species.out ./inputs/"${sample_input}"_reads.fastq.gz

Bracken-2.5/bracken -d db_kraken_"${db_input}"_species \
  -i "${sample_input}"_"${db_input}"_species.report \
  -o "${sample_input}"_"${db_input}"_species.bracken \
  --threads "$MAX_CPU"
