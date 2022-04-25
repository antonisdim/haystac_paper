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

mkdir -p db_krakenuniq_"${test_set}"

# allow kraken2 to use OMP multithreading
unset OMP_NUM_THREADS

# benchmark kraken building db performance test for 100 species
bin/krakenuniq-build --download-taxonomy --db db_krakenuniq_"${test_set}" \
  --threads "$MAX_CPU"

bin/krakenuniq-build \
  --add-to-library db_mutlifasta_inputs/db_input_"${test_set}".fasta \
  --db db_krakenuniq_"${test_set}" \
  --threads "$MAX_CPU"

bin/krakenuniq-build \
  --add-to-library krakenuniq_files/seqid2taxid.map \
  --db db_krakenuniq_"${test_set}" \
  --threads "$MAX_CPU"

bin/krakenuniq-build --build \
  --db db_krakenuniq_"${test_set}" \
  --taxids-for-genomes --taxids-for-sequences \
  --threads "$MAX_CPU" \
  --jellyfish-hash-size 30M

