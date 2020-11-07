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
use_conda=${args[2]}

#deactivating conda
haystack config \
  --use-conda "${use_conda}"

# run haystack sample analysis performance test for against a db of 5638 species without using conda

# delete any existing sample outputs so we can rebuild them
rm -rf ./samples/"${sample_input}"

# benchmark haystack sample analysis performance for samples of various sizes against a db of 5638 species

haystack sample \
  --sample-prefix "${sample_input}" \
  --fastq inputs/"${sample_input}".fastq.gz \
  --output ./samples/"${sample_input}"

haystack analyse \
  --mode abundances \
  --database ./haystack_db_"${db_input}"_input_8000_mem_conda_"${use_conda}" \
  --sample ./samples/"${sample_input}" \
  --output ./haystack_out_db_"${db_input}"_"${sample_input}"_conda_"${use_conda}"
