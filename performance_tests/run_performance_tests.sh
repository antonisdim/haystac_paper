#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# activate the conda environment for running the tests
eval "$(conda shell.bash hook)"
conda activate small_test

# bash strict mode (set after conda is activated, as the start-up scripts are not strict safe)
set -euo pipefail

# force bash to use the `time` command and not it's own native implementation
time=$(which time)

MAX_MEM=$(free -m | awk '/^Mem:/{printf "%.0f", $2}')

# list of database sizes (i.e. number of species) to test
species_counts=(10 100)

# make a folder to store the log files
mkdir -p logs

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for haystack                                                 -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes with a memory limit
for species_count in "${species_counts[@]}"; do

  for mem in $MAX_MEM; do
    echo "Running 'run_haystack_db_mem.sh' for ${species_count} species and ${mem} RAM"

    $time -v -o "logs/haystack_db_mem-${species_count}_species-${mem}_mem.time.log" \
      bash scripts/run_haystack_db_mem.sh "${species_count}_species" "${mem}" \
      >"logs/haystack_db_mem-${species_count}_species-${mem}_mem.log"
  done
done


# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_haystack_samples.sh' for input_10K reads, against ${species_count} species and use_conda True"

  $time -v -o "logs/haystack_samples-input_10K_reads-${species_count}_species-conda_True.time.log" \
    bash scripts/run_haystack_samples.sh input_10K_reads "${species_count}_species" "True" \
    >"logs/haystack_samples-input_10K_reads-${species_count}_species-conda_True.log"
done