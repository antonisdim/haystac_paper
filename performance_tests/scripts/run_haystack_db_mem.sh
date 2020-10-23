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

# delete any existing indices outputs so we can rebuild them
rm ../haystack_genome_cache/*/*.bt2l

# run haystack building db performance test with a mem limit for dbs of various sizes
haystack database \
  --mode build \
  --accessions ./haystack_configs/"haystack_db_${db_input}_input.txt" \
  --output ./"haystack_db_${db_input}_input_${mem_limit}_mem" \
  --mem "${mem_limit}" \
  --force-accessions
