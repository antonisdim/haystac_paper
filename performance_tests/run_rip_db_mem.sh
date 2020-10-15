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

haystack config \
  --email antonisdim41@gmail.com \
  --genome-cache-folder ../rip_genome_cache/

# delete any existing indices outputs so we can rebuild them
rm ../rip_genome_cache/*/*.bt2l

# run haystack building db performance test with a mem limit for dbs of various sizes
haystack database \
  -a ./haystack_configs/"rip_db_${test_set}_input.txt" \
  -o ./"rip_db_${test_set}_input_no_mem" \
  --mem $MEM_LIMIT_TEST

