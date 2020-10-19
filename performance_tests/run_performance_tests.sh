#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run all the tests for kraken/bracken

# setup a list of the different test_sets to run
test_sets_species=(10_species 100_species 1000_species 500_species 5638_species)
for test_set in "${test_sets_species[@]}"; do
  # building databases of different sizes
  /usr/bin/time -v bash run_kraken_db.sh "${test_set}" &>kraken_db.log
done

# setup a list of the different test_sets to run
test_sets_samples=(input_10K input_100K input_1M input_10M input_100M)
for test_set in "${test_sets_samples[@]}"; do
# analysing samples of different sizes against a db of 5638 species
  /usr/bin/time -v bash run_kraken_samples.sh "${test_set}" &>kraken_samples.log
done

# setup a list of the different test_sets to run
test_sets_number=(10 100 500 1000)
for test_set in "${test_sets_number[@]}"; do
  # analysing a sample of 1M reads against dbs of various sizes
  /usr/bin/time -v bash run_kraken_1M_vs_dbs.sh "${test_set}" &>kraken_1M_vs_dbs.log
done

# run all the tests for haystack

for test_set in "${test_sets_species[@]}"; do
  # building databases of different sizes with a memory limit
  /usr/bin/time -v bash run_haystack_db_mem.sh "${test_set}" &>rip_db_mem.log
done

# setup a list of the different test_sets to run
test_sets_samples_hsk=(input_10K_reads input_100K_reads input_1M_reads input_10M_reads input_100M_reads)
for test_set in "${test_sets_samples_hsk[@]}"; do
  # analysing samples of different sizes against a db of 5638 species with conda as a package manager
  /usr/bin/time -v bash run_haystack_samples_with_conda.sh "${test_set}" &>rip_samples_conda.log
done

# setup a list of the different test_sets to run
test_sets_1M_hsk=(10_species 100_species 500_species 1000_species)
for test_set in "${test_sets_1M_hsk[@]}"; do
  # analysing a sample of 1M reads against dbs of various sizes
  /usr/bin/time -v bash run_haystack_1M_vs_dbs.sh "${test_set}" &> rip_1M_vs_dbs.log
done

for test_set in "${test_sets_species[@]}"; do
  # building databases of different sizes with no memory limit
  /usr/bin/time -v bash run_haystack_db_no_mem.sh "${test_set}" &>rip_db_no_mem.log
done

for test_set in "${test_sets_samples_hsk[@]}"; do
  # analysing samples of different sizes against a db of 5638 species with no conda
  /usr/bin/time -v bash run_haystack_samples_no_conda.sh "${test_set}" &>rip_samples_no_conda.log
done


# run all the tests for sigma

for test_set in "${test_sets_species[@]}"; do
  # building databases of different sizes
  /usr/bin/time -v bash run_sigma_db.sh "${test_set}" &>sigma_db.log
done

for test_set in "${test_sets_samples[@]}"; do
  # analysing samples of different sizes against a db of 5638 species
  /usr/bin/time -v bash run_sigma_samples.sh "${test_set}" &>sigma_samples.log
done

# setup a list of the different test_sets to run
test_sets_1M_sigma=(1M_10sp 1M_100sp 1M_500sp 1M_1000sp)
for test_set in "${test_sets_1M_sigma[@]}"; do
  # analysing a sample of 1M reads against dbs of various sizes
  /usr/bin/time -v bash run_sigma_1M_vs_dbs.sh "${test_set}" &>sigma_1M_vs_dbs.log
done


# run all the tests for malt/hops

for test_set in "${test_sets_species[@]}"; do
  # building databases of different sizes
  /usr/bin/time -v bash run_malt_db.sh "${test_set}" &>malt_db.log
done

for test_set in "${test_sets_samples_hsk[@]}"; do
  # analysing samples of different sizes against a db of 5638 species
  /usr/bin/time -v bash run_malt_samples.sh "${test_set}" &>malt_samples.log
done

for test_set in "${test_sets_number[@]}"; do
  # analysing a sample of 1M reads against dbs of various sizes
  /usr/bin/time -v bash run_malt_1M_vs_dbs.sh "${test_set}" &>malt_1M_vs_dbs.log
done
