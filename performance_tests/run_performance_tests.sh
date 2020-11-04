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

# list of library sizes to test
read_counts=(input_10K input_100K)

# make a folder to store the log files
mkdir -p logs

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for haystack                                                 -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes with a memory limit
for species_count in "${species_counts[@]}"; do

  for mem in 8000 $MAX_MEM; do
    echo "Running 'run_haystack_db_mem.sh' for ${species_count} species and ${mem} RAM"

    # delete any existing indices outputs so we can rebuild them
    rm -f ./haystack_genome_cache_small/ncbi/*/*.bt2l
    rm -f ./haystack_genome_cache_small/ncbi/*/*.fasta.gz.fai

    $time -v -o "logs/haystack_db_mem-${species_count}_species-${mem}_mem.time.log" \
      bash scripts/run_haystack_db_mem.sh "${species_count}_species" "${mem}" \
      >"logs/haystack_db_mem-${species_count}_species-${mem}_mem.log"
  done
done

# analysing samples of different sizes against a db of 5636 species with conda as a package manager
for read_count in "${read_counts[@]}"; do

  for bool in "False" "True"; do
    echo "Running 'run_haystack_samples.sh' for ${read_count} reads, against 5636 species and use_conda ${bool}"

    $time -v -o "logs/haystack_samples-${read_count}_reads-5636_species-conda_${bool}.time.log" \
      bash scripts/run_haystack_samples.sh "${read_count}_reads" 5636_species "${bool}" \
      >"logs/haystack_samples-${read_count}_reads-5636_species-conda_${bool}.log"
  done
done
