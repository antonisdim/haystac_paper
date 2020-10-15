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

# benchmark malt sample analysis performance test for a sample of 1M reads against dbs of various size
malt-run \
  -i ./inputs/input_1M_reads.fastq \
  -d index_new_"${test_set}"_species \
  -m BlastN \
  -o input_1M_reads_"${test_set}"sp \
  -t "$MAX_CPU"
hops --mode me_po \
  --input input_1M_reads_"${test_set}"sp.rma6 \
  --output hops_input_1M_reads_"${test_set}"sp \
  --configFile configfile.txt
