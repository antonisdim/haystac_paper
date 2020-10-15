#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# create a new conda environment to run all the test_sets
eval "$(conda shell.bash hook)"
conda env create -f environment.yaml -n performance_test
conda activate performance_test

# download the binaries for Bracken and Sigma
wget https://github.com/jenniferlu717/Bracken/archive/v2.5.tar.gz -O Bracken-2.5.tar.gz
tar xvzf Bracken-2.5.tar.gz
rm Bracken-2.5.tar.gz

wget --content-disposition http://sourceforge.net/projects/sigma-omicsbio/files/V1.0.3/sigma-v1.0.3.tar.gz/download \
  -O Sigma.tar.gz
tar xvzf Sigma.tar.gz
rm Sigma.tar.gz

# setup a list of the different test_sets to run
test_sets=(10_species 100_species 500_species 1000_species 5638_species)

# get the max available memory (in GB)
MAX_MEM=$(free -m | awk '/^Mem:/{printf "%.0f", $2/1024}')

# configure MALT to use the max available memory
sed -i'' "s/-Xmx.*/-Xmx${MAX_MEM}G/" ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-run.vmoptions
sed -i'' "s/-Xmx.*/-Xmx${MAX_MEM}G/" ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-build.vmoptions

# create mutlifasta files for kraken2 and MALT and move them to the right folders

mkdir -p db_mutlifasta_inputs

for test_set in "${test_sets[@]}"; do
  cat "${test_set}.txt" | xargs zcat >"db_mutlifasta_inputs/db_input_${test_set}.fasta"
done

# download mapping for MALT
wget -P mapping_files https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Jul2020-2.db.zip
wget -P mapping_files https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jul2020.db.zip
gunzip -c mapping_files/megan-map-Jul2020-2.db.zip >mapping_files/megan-map-Jul2020-2.db

# make the database folders for Sigma by symlinking the genomes into the relevant folders
for test_set in "${test_sets[@]}"; do
  mkdir -p "sigma_db_${test_set}"
  cat "${test_set}.txt" | xargs dirname | awk '{print "../"$0}' | xargs ln -s -t "sigma_db_${test_set}"
done

# run the performance test_sets
bash run_performance_test.sh &>performance_test.log
