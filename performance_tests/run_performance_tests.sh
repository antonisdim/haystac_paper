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

MAX_MEM=$(free -m | awk '/^Mem:/{printf "%.0f", $2/1024}')

# list of database sizes (i.e. number of species) to test
species_counts=(10 100 500 1000 5638)

# list of library sizes to test
read_counts=(input_10K input_100K input_1M input_10M input_100M)

# make a folder to store the log files
mkdir -p logs

# ----------------------------------------------------------------------------------------------------------------------
# -                                   run all the tests for kraken/bracken                                             -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
  echo "Running 'run_kraken_db.sh' for ${species_count} species"

  /usr/bin/time -v -o "logs/kraken_db_${species_count}_species.time.log" \
    bash run_kraken_db.sh "${species_count}_species" &>"logs/kraken_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_kraken_samples.sh' for ${read_count} reads against 5638 species"

  /usr/bin/time -v -o "logs/kraken_samples_${read_count}.time.log" \
    bash run_kraken_samples.sh "${read_count}" 5638 &>"logs/kraken_samples_${read_count}.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_kraken_samples.sh' for input_1M reads against ${species_count} species"

  /usr/bin/time -v -o "logs/kraken_1M_vs_dbs_${species_count}_species.time.log" \
    bash run_kraken_samples.sh input_1M "${species_count}" &>"logs/kraken_1M_vs_dbs_${species_count}.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for haystack                                                 -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes with a memory limit
for species_count in "${species_counts[@]}"; do

  for mem in 8000 $MAX_MEM; do
    echo "Running 'run_haystack_db_mem.sh' for ${species_count} species and ${mem} RAM"

    /usr/bin/time -v -o "logs/rip_db_mem_${species_count}_species_${mem}_mem.time.log" \
      bash run_haystack_db_mem.sh "${species_count}_species" "${mem}" &>"logs/rip_db_${mem}_mem_${species_count}.log"
  done
done

# analysing samples of different sizes against a db of 5638 species with conda as a package manager
for read_count in "${read_counts[@]}"; do

  for bool in "False" "True"; do
    echo "Running 'run_haystack_samples.sh' for ${read_count} reads, against 5638 species and use_conda ${bool}"

    /usr/bin/time -v -o "logs/rip_samples_conda_${read_count}_reads_5638_species_conda_${bool}.time.log" \
      bash run_haystack_samples.sh "${read_count}_reads" 5638_species "${bool}" \
      &>"logs/rip_samples_conda_${read_count}_5638_species_conda_${bool}.log"
  done
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_haystack_samples.sh' for input_1M reads, against ${species_count} species and use_conda True"

  /usr/bin/time -v -o "logs/rip_1M_vs_dbs_${species_count}_species.time.log" \
    bash run_haystack_samples.sh input_1M_reads "${species_count}_species" "True" \
    &>"logs/rip_1M_vs_dbs_${species_count}_conda_True.log"
done




# ----------------------------------------------------------------------------------------------------------------------
# -                                       run all the tests for sigma                                                  -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
   echo "Running 'run_sigma_db.sh' for ${species_count} species"

  /usr/bin/time -v -o "logs/sigma_db_${species_count}_species.time.log" \
    bash run_sigma_db.sh "${species_count}_species" &>"logs/sigma_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_sigma_samples.sh' for ${read_count} reads against 5638 species"

  /usr/bin/time -v -o "logs/sigma_samples_${read_count}_5638_species.time.log" \
    bash run_sigma_samples.sh "${read_count}" 5638 &>"logs/sigma_samples_${read_count}_reads_5638_species.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_sigma_samples.sh' for input_1M reads against ${species_count} species"

  /usr/bin/time -v -o "logs/sigma_1M_vs_dbs_${species_count}_species.time.log" \
    bash run_sigma_samples.sh input_1M "${species_count}" &>"logs/sigma_1M_vs_dbs_${species_count}.log"
done

# ----------------------------------------------------------------------------------------------------------------------
# -                                     run all the tests for malt/hops                                                -
# ----------------------------------------------------------------------------------------------------------------------

# building databases of different sizes
for species_count in "${species_counts[@]}"; do
  echo "Running 'run_malt_db.sh' for ${species_count} species"

  /usr/bin/time -v -o "logs/malt_db_${species_count}_species.time.log" \
    bash run_malt_db.sh "${species_count}_species" &>"logs/malt_db_${species_count}.log"
done

# analysing samples of different sizes against a db of 5638 species
for read_count in "${read_counts[@]}"; do
  echo "Running 'run_malt_samples.sh' for ${read_count} reads against 5638 species"

  /usr/bin/time -v -o "logs/malt_samples_${read_count}_reads_5638_species.time.log" \
    bash run_malt_samples.sh "${read_count}_reads" 5638 &>"logs/malt_samples_${read_count}_5638_species.log"
done

# analysing a sample of 1M reads against dbs of various sizes
for species_count in "${species_counts[@]:0:4}"; do
  echo "Running 'run_malt_samples.sh' for input_1M reads against ${species_count} species"

  /usr/bin/time -v -o "logs/malt_1M_vs_dbs_1${species_count}_species.time.log" \
    bash run_malt_samples.sh input_1M_reads "${species_count}" &>"logs/malt_1M_vs_dbs_${species_count}.log"
done
