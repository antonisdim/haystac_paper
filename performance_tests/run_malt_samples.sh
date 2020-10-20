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

# benchmark malt sample analysis performance test for samples of 1M reads against dbs of various size
malt-run \
  -t $MAX_CPU \
  -i ./inputs/"${test_set}".fastq \
  -d index_new_5638_species \
  -m BlastN \
  -o .
hops --mode me_po \
  --input "${test_set}".rma6 \
  --output hops_"${test_set}" \
  --configFile configfile.txt
