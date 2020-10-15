#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

args=("$@")
test_set=${args[0]}


# delete any existing indices outputs so we can rebuild them
rm ../rip_genome_cache/*/*.bt2l

# run haystack building db performance test with no mem limit for dbs of various sizes
haystack database \
  -a ./haystack_configs/"rip_db_${test_set}_input.txt" \
  -o ./"rip_db_${test_set}_input_no_mem"

