#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# run sigma building db perfromance test of dbs of various sizes

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

for i in ./sigma_db_10_species/*; do rm $i/*.bt2; done
time -v ./Sigma/bin/sigma-index-genomes -c sigma_config_10_species.cfg -w ./sigma_db_10_species -p $MAX_CPU
for i in ./sigma_db_100_species/*; do rm $i/*.bt2; done
time -v ./Sigma/bin/sigma-index-genomes -c sigma_config_100_species.cfg -w ./sigma_db_100_species -p $MAX_CPU
for i in ./sigma_db_1000_species/*; do rm $i/*.bt2; done
time -v ./Sigma/bin/sigma-index-genomes -c sigma_config_1000_species.cfg -w ./sigma_db_1000_species -p $MAX_CPU
for i in ./sigma_db_500_species/*; do rm $i/*.bt2; done
time -v ./Sigma/bin/sigma-index-genomes -c sigma_config_500_species.cfg -w ./sigma_db_500_species -p $MAX_CPU
for i in ./sigma_db_5638_species/*; do rm $i/*.bt2; done
time -v ./Sigma/bin/sigma-index-genomes -c sigma_config_5638_species.cfg -w ./sigma_db_5638_species -p $MAX_CPU
