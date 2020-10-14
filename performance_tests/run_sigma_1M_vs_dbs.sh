#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail


MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

mkdir sigma_outputs/input_1M_10sp
mkdir sigma_outputs/input_1M_100sp
mkdir sigma_outputs/input_1M_1000sp
mkdir sigma_outputs/input_1M_500sp

#run sigma sample analysis performance test for a sample of 1M reads against a db of 10 species
/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c ./sigma_configs/sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c ./sigma_configs/sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp/; \
	./Sigma/bin/sigma-solve-model \
	-c ./sigma_configs/sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_10sp/

#run sigma sample analysis performance test for a sample of 1M reads against a db of 100 species
/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c ./sigma_configs/sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c ./sigma_configs/sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp/; \
	./Sigma/bin/sigma-solve-model \
	-c ./sigma_configs/sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_100sp/

#run sigma sample analysis performance test for a sample of 1M reads against a db of 1000 species
/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c ./sigma_configs/sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c ./sigma_configs/sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp/; \
	./Sigma/bin/sigma-solve-model \
	-c ./sigma_configs/sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_1000sp/

#run sigma sample analysis performance test for a sample of 1M reads against a db of 500 species
/usr/bin/time -v ./Sigma/bin/sigma-align-reads \
	-c ./sigma_configs/sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c ./sigma_configs/sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp/; \
	./Sigma/bin/sigma-solve-model \
	-c ./sigma_configs/sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_500sp/

