#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# activate the conda environment for running the tests
eval "$(conda shell.bash hook)"
conda activate performance_test

# TODO move all other other scripts into a script folder (only leave setup_tests and run_tests in this folder)

# TODO please consolidate all these lists. There should only be one instances of each list.

# list of database sizes (i.e. number of species) to test
species_counts=(10 100 500 1000)

# TODO no need to fixed string suffixes here, move `_species` into the command
species_counts2=(10_species 100_species 500_species 1000_species)

# TODO just append `5638` onto the earlier list
species_counts_all=(10_species 100_species 500_species 1000_species 5638_species)

# list of library sizes to test
read_counts=(input_10K input_100K input_1M input_10M input_100M)

# TODO refactor this out...
read_counts2=(input_10K_reads input_100K_reads input_1M_reads input_10M_reads input_100M_reads)

# TODO why are you using so many redundant lists?
species_counts_sigma=(1M_10sp 1M_100sp 1M_500sp 1M_1000sp)

# make a folder to store the log files
mkdir -p logs

# ----------------------------------------------------------------------------------------------------------------------
# -                                   run all the tests for kraken/bracken                                             -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts_all[@]}"; do
  # TODO add progress messages to every for loop, so we know what's going on when we run the script, e.g.
  echo "Running 'run_kraken_db.sh' for ${species_count} species"

  /usr/bin/time -v -o "logs/kraken_db_${species_count}.time.log" \
    bash run_kraken_db.sh "${species_count}" &>"logs/kraken_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts[@]}"; do
  /usr/bin/time -v -o "logs/kraken_samples_${read_count}.time.log" \
    bash run_kraken_samples.sh "${read_count}" &>"logs/kraken_samples_${read_count}.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]}"; do
  /usr/bin/time -v -o "logs/kraken_1M_vs_dbs_${species_count}.time.log" \
    bash run_kraken_1M_vs_dbs.sh "${species_count}" &>"logs/kraken_1M_vs_dbs_${species_count}.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for haystack                                                 -
# ----------------------------------------------------------------------------------------------------------------------

# TODO we need to run all the `haystack` tests in
#   1. restricted memory mode, w/ conda
#   2. restricted memory mode, w/ no-conda
#   3. unrestricted memory, w/ conda
#   4. unrestricted memory, w/ no-conda
#
#   make two extra loops, e.g. `for mem_limit in 8000, MAX_MEM` and pass params into each script
#   there should only be one `run_haystack.sh` script! same goes for ALL the other methods

# building databases of different sizes with a memory limit
for species_count in "${species_counts_all[@]}"; do
  /usr/bin/time -v -o "logs/rip_db_mem_${species_count}.time.log" \
    bash run_haystack_db_mem.sh "${species_count}" &>"logs/rip_db_mem_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species with conda as a package manager
for read_count in "${read_counts2[@]}"; do
  /usr/bin/time -v -o "logs/rip_samples_conda_${read_count}.time.log" \
    bash run_haystack_samples_with_conda.sh "${read_count}" &>"logs/rip_samples_conda_${read_count}.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_counts in "${species_counts2[@]}"; do
  /usr/bin/time -v -o "logs/rip_1M_vs_dbs_${species_count}.time.log" \
    bash run_haystack_1M_vs_dbs.sh "${species_counts}" &>"logs/rip_1M_vs_dbs_${species_counts}.log"
done

# building databases of different sizes with no memory limit
for species_count in "${species_counts_all[@]}"; do
  /usr/bin/time -v -o "logs/rip_db_no_mem_${species_count}.time.log" \
    bash run_haystack_db_no_mem.sh "${species_count}" &>"logs/rip_db_no_mem_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species with no conda
for read_count in "${read_counts2[@]}"; do
  /usr/bin/time -v -o "logs/rip_samples_no_conda_${read_count}.time.log" \
    bash run_haystack_samples_no_conda.sh "${read_count}" &>"logs/rip_samples_no_conda_${read_count}.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                       run all the tests for sigma                                                  -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts_all[@]}"; do
  /usr/bin/time -v -o "logs/sigma_db_${read_count}.time.log" \
    bash run_sigma_db.sh "${species_count}" &>"logs/sigma_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts[@]}"; do
  /usr/bin/time -v -o "logs/sigma_samples_${read_count}.time.log" \
    bash run_sigma_samples.sh "${read_count}" &>"logs/sigma_samples_${read_count}.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts_sigma[@]}"; do
  /usr/bin/time -v -o "logs/sigma_1M_vs_dbs_${species_count}.time.log" \
    bash run_sigma_1M_vs_dbs.sh "${species_count}" &>"logs/sigma_1M_vs_dbs_${species_count}.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for malt/hops                                                -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts_all[@]}"; do
  /usr/bin/time -v -o "logs/malt_db_${species_count}.time.log" \
    bash run_malt_db.sh "${species_count}" &>"logs/malt_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts2[@]}"; do
  /usr/bin/time -v -o "logs/malt_samples_${read_count}.time.log" \
    bash run_malt_samples.sh "${read_count}" &>"logs/malt_samples_${read_count}.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]}"; do
  /usr/bin/time -v -o "logs/malt_1M_vs_dbs_${species_count}.time.log" \
    bash run_malt_1M_vs_dbs.sh "${species_count}" &>"logs/malt_1M_vs_dbs_${species_count}.log"
done
