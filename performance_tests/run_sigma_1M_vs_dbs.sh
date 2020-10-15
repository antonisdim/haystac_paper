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

mkdir -p sigma_outputs/"input_${test_set}"

# benchmark Sigma sample analysis performance test for a sample of 1M reads for dbs of various sizes
./Sigma/bin/sigma-align-reads \
  -c ./sigma_configs/sigma_config_input_"${test_set}".cfg \
  -w sigma_outputs/input_"${test_set}" \
  -p "$MAX_CPU"
./Sigma/bin/sigma-build-model \
  -c ./sigma_configs/sigma_config_input_"${test_set}".cfg \
  -w sigma_outputs/input_"${test_set}"/
./Sigma/bin/sigma-solve-model \
  -c ./sigma_configs/sigma_config_input_"${test_set}".cfg \
  -w sigma_outputs/input_"${test_set}"/ \
  -t "$MAX_CPU"
mv sigma_out.* ./sigma_outputs/input_"${test_set}"/


