#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# fetch genomes for the performance tests and time it

haystack config \
  --email antonisdim41@gmail.com \
  --cache ./rip_genome_cache/

haystack database --snakemake '{"conda_create_envs_only": 1}'
haystack sample --snakemake '{"conda_create_envs_only": 1}'
haystack analyse --snakemake '{"conda_create_envs_only": 1}'

time -v haystack database \
  --mode fetch \
  --accessions ./haystack_configs/"rip_db_5638_species_input.txt" \
  --output ./"genome_fetch"
