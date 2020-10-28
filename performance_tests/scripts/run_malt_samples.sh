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

# benchmark malt sample analysis performance test for samples of 1M reads against dbs of various size
malt-run \
  -t $MAX_CPU \
  -i ./inputs/"${sample_input}".fastq \
  -d index_new_"${db_input}"_species \
  -m BlastN \
  -o "${sample_input}"_"${db_input}"_species.rma6

hops --mode me_po \
  --input "${sample_input}"_"${db_input}"_species.rma6 \
  --output hops_"${sample_input}"_"${db_input}"_species \
  --configFile hops_files/configfile.txt


