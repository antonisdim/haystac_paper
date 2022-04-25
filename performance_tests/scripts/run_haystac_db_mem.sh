#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
db_input=${args[0]}
mem_limit=${args[1]}
use_conda=${args[2]}

#deactivating conda

haystac config \
  --use-conda "${use_conda}"

# run haystac building db performance test with a mem limit for dbs of various sizes

haystac database \
  --mode build \
  --accessions ./haystac_configs/"haystac_db_${db_input}_refseq_input.txt" \
  --output ./"haystac_db_${db_input}_input_${mem_limit}_mem_conda_${use_conda}" \
  --mem "${mem_limit}" \
  --force-accessions
