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
bin/krakenuniq \
  --db db_krakenuniq_"${db_input}"_species \
  --threads "$MAX_CPU" \
  --fastq-input \
  --gzip-compressed \
  --report-file "${sample_input}"_"${db_input}"_species_krakenuniq.report \
  --unclassified-out "${sample_input}"_"${db_input}"_species_unclassified_krakenuniq.out \
  --classified-out "${sample_input}"_"${db_input}"_species_classified_krakenuniq.out \
  --output "${sample_input}"_"${db_input}"_species_krakenuniq.out ./inputs/"${sample_input}"_reads.fastq.gz
