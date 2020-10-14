#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run haystack sample analysis performance test for a sample of 1M reads against a db of 10 species

# setup a list of the different test_sets to run
test_sets=(10_species 100_species 500_species 1000_species)

for test_set in "${test_sets[@]}"; do
  # delete any existing sample outputs so we can rebuild them
  rm -r ./samples/input_1M_reads

  # run haystack sample analysis performance test for a sample of 1M reads for dbs of various sizes
  /usr/bin/time -v haystack analyse \
	  --mode abundances \
	  --database ./"rip_db_${test_set}_input" \
	  --sample ./samples/input_1M_reads \
	  --output ./"rip_db_${test_set}_analysis_output"

done
