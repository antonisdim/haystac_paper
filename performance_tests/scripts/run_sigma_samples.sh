#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
sample_input=${args[0]}
db_input=${args[1]}

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

mkdir -p sigma_outputs

mkdir -p sigma_outputs/"${sample_input}"

# benchmark Sigma sample analysis performance for samples of various sizes against a db of 5638 species
./Sigma/bin/sigma-align-reads \
  -c ./sigma_configs/sigma_config_"${sample_input}"_"${db_input}"_sp.cfg \
  -w sigma_outputs/"${sample_input}" \
  -p "$MAX_CPU"

./Sigma/bin/sigma-build-model \
  -c ./sigma_configs/sigma_config_"${sample_input}"_"${db_input}"_sp.cfg \
  -w ./sigma_outputs/"${sample_input}"/

./Sigma/bin/sigma-solve-model \
  -t "$MAX_CPU" \
  -c ./sigma_configs/sigma_config_"${sample_input}"_"${db_input}"_sp.cfg \
  -w ./sigma_outputs/"${sample_input}"

mv sigma_out.* ./sigma_outputs/"${sample_input}"/
