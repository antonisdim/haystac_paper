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

# TODO remove the sigma and bracken binaries from the git repo and download them here instead
# TODO this can be in the environment.yaml right?
# install haystack
python -m pip install git+https://github.com/antonisdim/haystack

# setup a list of the different test_sets to run
test_sets=(10_species 100_species 500_species 1000_species 5638_species)

# get the max available memory (in GB)
MAX_MEM=$(free -m | awk '/^Mem:/{printf "%.0f", $2/1024}')

# configure MALT to use the max available memory
sed -i'' "s/-Xmx.*/-Xmx${MAX_MEM}G/" ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-run.vmoptions
sed -i'' "s/-Xmx.*/-Xmx${MAX_MEM}G/" ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-build.vmoptions

# create mutlifasta files for kraken2 and MALT and move them to the right folders
for test_set in "${test_sets[@]}"; do
  xargs <"${test_set}.txt" | cat >"db_input_${test_set}.fasta.gz"
done

# create a mutlifasta file for all genomes in the database
cat ../rip_genome_cache/*/*fasta.gz >db_input_5638_species.fasta.gz

# TODO do we really need both *.fasta and *.fasta.gz?
gunzip -d db_input_1000_species.fasta.gz
gunzip -d db_input_100_species.fasta.gz
gunzip -d db_input_10_species.fasta.gz
gunzip -d db_input_500_species.fasta.gz
gunzip -d db_input_5638_species.fasta.gz

mkdir db_mutlifasta_inputs
mv ./*.fasta db_mutlifasta_inputs/

# download mapping for MALT
wget -P mapping_files https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Jul2020-2.db.zip
wget -P mapping_files https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jul2020.db.zip
gunzip -c mapping_files/megan-map-Jul2020-2.db.zip >mapping_files/megan-map-Jul2020-2.db

# make the database folders for Sigma by symlinking the genomes into the relevant folders
for test_set in "${test_sets[@]}"; do
  mkdir -p "sigma_db_${test_set}"
  xargs dirname <"${test_set}.txt" | xargs ln -s -t "sigma_db_${test_set}"
done

# run the performance test_sets
bash run_performance_test.sh &>performance_test.log
