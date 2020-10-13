#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run sigma sample analysis performance test for samples of various size

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

mkdir sigma_outputs/input_10K
mkdir sigma_outputs/input_100K
mkdir sigma_outputs/input_1M
mkdir sigma_outputs/input_10M
mkdir sigma_outputs/input_100M

/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_10K.cfg \
	-w sigma_outputs/input_10K \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_10K.cfg \
	-w ./sigma_outputs/input_10K/; \
	./Sigma/bin/sigma-solve-model \
	-t $MAX_CPU \
	-c sigma_config_input_10K \
	-w ./sigma_outputs/input_10K/; \
	mv sigma_out.* ./sigma_outputs/input_10K/

/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_100K.cfg \
	-w sigma_outputs/input_100K \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_100K.cfg \
	-w ./sigma_outputs/input_100K/; \
	./Sigma/bin/sigma-solve-model \
	-t $MAX_CPU \
	-c sigma_config_input_100K \
	-w ./sigma_outputs/input_100K/; \
	mv sigma_out.* ./sigma_outputs/input_100K/

/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_1M.cfg \
	-w sigma_outputs/input_1M \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_1M.cfg \
	-w ./sigma_outputs/input_1M/; \
	./Sigma/bin/sigma-solve-model \
	-t $MAX_CPU \
	-c sigma_config_input_1M \
	-w ./sigma_outputs/input_1M/; \
	mv sigma_out.* ./sigma_outputs/input_1M/

/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_10M.cfg \
	-w sigma_outputs/input_10M \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_10M.cfg \
	-w ./sigma_outputs/input_10M/; \
	./Sigma/bin/sigma-solve-model \
	-t $MAX_CPU \
	-c sigma_config_input_10M \
	-w ./sigma_outputs/input_10M/; \
	mv sigma_out.* ./sigma_outputs/input_10M/

/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_100M.cfg \
	-w sigma_outputs/input_100M \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_100M.cfg \
	-w ./sigma_outputs/input_100M/; \
	./Sigma/bin/sigma-solve-model \
	-t $MAX_CPU \
	-c sigma_config_input_100M \
	-w ./sigma_outputs/input_100M/; \
	mv sigma_out.* ./sigma_outputs/input_100M/

