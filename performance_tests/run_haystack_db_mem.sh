#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
test_set=${args[0]}

MEM_LIMIT_TEST=8000

# delete any existing indices outputs so we can rebuild them
rm ../rip_genome_cache/*/*.bt2l

# run haystack building db performance test with a mem limit for dbs of various sizes
haystack database \
  --mode build \
  --accessions ./haystack_configs/"haystack_db_${test_set}_input.txt" \
  --output ./"haystack_db_${test_set}_input_mem" \
  --mem $MEM_LIMIT_TEST

