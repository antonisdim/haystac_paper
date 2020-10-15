#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run sigma building db performance test of dbs of various sizes

args=("$@")
test_set=${args[0]}

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# delete any existing indices so we can rebuild them
rm ./sigma_db_"${test_set}"/*/*.bt2

# benchmark Sigma
./Sigma/bin/sigma-index-genomes \
  -c ./sigma_configs/"sigma_config_${test_set}.cfg" \
  -w "./sigma_db_${test_set}" \
  -p "${MAX_CPU}"

