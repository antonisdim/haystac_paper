#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# activate the conda environment for running the tests
eval "$(conda shell.bash hook)"
conda activate performance_test

# bash strict mode (set after conda is activated, as the start-up scripts are not strict safe)
set -euo pipefail

# force bash to use the `time` command and not it's own native implementation
time=$(which time)

MAX_MEM=$(free -m | awk '/^Mem:/{printf "%.0f", $2/1024}')

# list of database sizes (i.e. number of species) to test
species_counts=(10 100 500 1000 5616)

# list of library sizes to test
read_counts=(input_10K input_100K input_1M input_10M input_100M)

# make a folder to store the log files
mkdir -p logs

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for haystack                                                 -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes with a memory limit
for species_count in "${species_counts[@]}"; do

  for mem in 8000 $MAX_MEM; do
    echo "Running 'run_haystack_db_mem.sh' for ${species_count} species and ${mem} RAM"

    $time -v -o "logs/haystack_db_mem-${species_count}_species-${mem}_mem.time.log" \
      bash scripts/run_haystack_db_mem.sh "${species_count}_species" "${mem}" \
      &>"logs/haystack_db_mem-${species_count}_species-${mem}_mem.log"
  done
done

# analysing samples of different sizes against a db of 5616 species with conda as a package manager
for read_count in "${read_counts[@]}"; do

  for bool in "False" "True"; do
    echo "Running 'run_haystack_samples.sh' for ${read_count} reads, against 5616 species and use_conda ${bool}"

    $time -v -o "logs/haystack_samples-${read_count}_reads-5616_species-conda_${bool}.time.log" \
      bash scripts/run_haystack_samples.sh "${read_count}_reads" 5616_species "${bool}" \
      &>"logs/haystack_samples-${read_count}_reads-5616_species-conda_${bool}.log"
  done
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_haystack_samples.sh' for input_1M reads, against ${species_count} species and use_conda True"

  $time -v -o "logs/haystack_samples-input_1M_reads-${species_count}_species-conda_True.time.log" \
    bash scripts/run_haystack_samples.sh input_1M_reads "${species_count}_species" "True" \
    &>"logs/haystack_samples-input_1M_reads-${species_count}_species-conda_True.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                   run all the tests for kraken/bracken                                             -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
  echo "Running 'run_kraken_db.sh' for ${species_count} species"

  $time -v -o "logs/kraken_db-${species_count}_species.time.log" \
    bash scripts/run_kraken_db.sh "${species_count}_species" &>"kraken_db-${species_count}_species.log"
done

# analysing samples of different sizes against a db of 5616 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_kraken_samples.sh' for ${read_count} reads against 5616 species"

  $time -v -o "logs/kraken_samples-${read_count}_reads-5616_species.time.log" \
    bash scripts/run_kraken_samples.sh "${read_count}" 5616 \
    &>"logs/kraken_samples-${read_count}_reads-5616_species.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_kraken_samples.sh' for input_1M reads against ${species_count} species"

  $time -v -o "logs/kraken_samples-input_1M-${species_count}_species.time.log" \
    bash scripts/run_kraken_samples.sh input_1M "${species_count}" \
    &>"logs/kraken_samples-input_1M-${species_count}_species.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                       run all the tests for sigma                                                  -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
  echo "Running 'run_sigma_db.sh' for ${species_count} species"

  $time -v -o "logs/sigma_db-${species_count}_species.time.log" \
    bash scripts/run_sigma_db.sh "${species_count}_species" &>"logs/sigma_db-${species_count}_species.log"
done

# analysing samples of different sizes against a db of 5616 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_sigma_samples.sh' for ${read_count} reads against 5616 species"

  $time -v -o "logs/sigma_samples-${read_count}_reads-5616_species.time.log" \
    bash scripts/run_sigma_samples.sh "${read_count}" 5616 &>"logs/sigma_samples-${read_count}_reads-5616_species.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_sigma_samples.sh' for input_1M reads against ${species_count} species"

  $time -v -o "logs/sigma_samples-input_1M_reads-${species_count}_species.time.log" \
    bash scripts/run_sigma_samples.sh input_1M "${species_count}" &>"logs/sigma_samples-input_1M-${species_count}_species.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for malt/hops                                                -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
  echo "Running 'run_malt_db.sh' for ${species_count} species"

  $time -v -o "logs/malt_db-${species_count}_species.time.log" \
    bash scripts/run_malt_db.sh "${species_count}_species" &>"logs/malt_db-${species_count}_species.log"
done

# analysing samples of different sizes against a db of 5616 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_malt_samples.sh' for ${read_count} reads against 5616 species"

  $time -v -o "logs/malt_samples-${read_count}_reads-5616_species.time.log" \
    bash scripts/run_malt_samples.sh "${read_count}_reads" 5616 &>"logs/malt_samples-${read_count}_reads-5616_species.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_malt_samples.sh' for input_1M reads against ${species_count} species"

  $time -v -o "logs/malt_samples-input_1M_reads-${species_count}_species.time.log" \
    bash scripts/run_malt_samples.sh input_1M_reads "${species_count}" \
    &>"logs/malt_samples-input_1M_reads-${species_count}_species.log"
done
