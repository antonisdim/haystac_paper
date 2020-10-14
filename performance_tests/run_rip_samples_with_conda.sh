#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run haystack sample analysis performance test for against a db of 5638 species using conda

# setup a list of the different test_sets to run
test_sets=(input_10K_reads input_100K_reads input_100M_reads input_1M_reads input_10M_reads)

for test_set in "${test_sets[@]}"; do
  # delete any existing sample outputs so we can rebuild them
  rm -r ./samples/"${test_set}"

  # benchmark haystack sample analysis performance for samples of various sizes against a db of 5638 species

  /usr/bin/time -v haystack sample \
    --sample-prefix "${test_set}" \
    --fastq ./inputs/"${test_set}".fastq.gz \
    --output ./samples/"${test_set}"

  /usr/bin/time -v haystack analyse \
    --mode abundances \
    --database ./refseq_resp \
    --sample ./samples/"${test_set}" \
    --output ./refseq_analysis_output

done
