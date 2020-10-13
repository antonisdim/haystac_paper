#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run sigma building db performance test of dbs of various sizes

# setup a list of the different test_sets to run
test_sets=(10_species 100_species 500_species 1000_species 5638_species)

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

for test_set in "${test_sets[@]}"; do
  # delete any existing indices so we can rebuild them
  rm ./sigma_db_"${test_set}"/*/*.bt2

  # benchmark Sigma
  /usr/bin/time -v ./Sigma/bin/sigma-index-genomes -c "sigma_config_${test_set}.cfg" -w "./sigma_db_${test_set}" -p "${MAX_CPU}"
done
