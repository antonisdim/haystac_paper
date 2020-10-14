#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

mkdir -p sigma_outputs

# setup a list of the different test_sets to run
test_sets=(input_10K input_100K input_100M input_1M input_10M)

for test_set in "${test_sets[@]}"; do

  mkdir -p sigma_outputs/"${test_set}"

  # benchmark Sigma sample analysis performance for samples of various sizes against a db of 5638 species
  /usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	  -c ./sigma_configs/sigma_config_"${test_set}".cfg \
	  -w sigma_outputs/"${test_set}" \
	  -p "$MAX_CPU"; \
	  ./Sigma/bin/sigma-build-model \
	  -c ./sigma_configs/sigma_config_"${test_set}".cfg \
	  -w ./sigma_outputs/"${test_set}"/; \
	  ./Sigma/bin/sigma-solve-model \
	  -t "$MAX_CPU" \
	  -c ./sigma_configs/sigma_config_"${test_set}" \
	  -w ./sigma_outputs/"${test_set}"/; \
	  mv sigma_out.* ./sigma_outputs/"${test_set}"/


done

