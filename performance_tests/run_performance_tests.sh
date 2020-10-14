#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run all the tests for kraken/bracken

# building databases of different sizes
bash run_kraken_db.sh &> kraken_db.log
# analysing samples of different sizes against a db of 5638 species
bash run_kraken_samples.sh &> kraken_samples.log
# analysing a sample of 1M reads against dbs of various sizes
bash run_kraken_1M_vs_dbs.sh &> kraken_1M_vs_dbs.log

# run all the tests for haystack

# building databases of different sizes with a memory limit
bash run_rip_db_mem.sh &> rip_db_mem.log
# analysing samples of different sizes against a db of 5638 species with conda as a package manager
bash run_rip_samples.sh &> rip_samples.log
# analysing a sample of 1M reads against dbs of various sizes
bash run_rip_1M_vs_dbs.sh &> rip_1M_vs_dbs.log
# building databases of different sizes with no memory limit
bash run_rip_db_no_mem.sh &> rip_db_no_mem.log
# analysing samples of different sizes against a db of 5638 species with no conda
bash run_rip_samples_no_conda.sh &> rip_samples_no_conda.log

# run all the tests for sigma

# building databases of different sizes
bash run_sigma_db.sh &> sigma_db.log
# analysing samples of different sizes against a db of 5638 species
bash run_sigma_samples.sh &> sigma_samples.log
# analysing a sample of 1M reads against dbs of various sizes
bash run_sigma_1M_vs_dbs.sh &> sigma_1M_vs_dbs.log

# run all the tests for malt/hops

# building databases of different sizes
bash run_malt_db.sh &> malt_db.log
# analysing samples of different sizes against a db of 5638 species
bash run_malt_samples.sh &> malt_samples.log
# analysing a sample of 1M reads against dbs of various sizes
bash run_malt_1M_vs_dbs.sh &> malt_1M_vs_dbs.log
