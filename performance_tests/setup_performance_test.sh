#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

eval "$(conda shell.bash hook)"
# bash strict mode
set -euo pipefail

conda env create -f environment.yaml -n performance_test

# activate conda env 

conda activate performance_test

# set the right options for ram usage for malt 

sed -i'' 's/-Xmx64G/-Xmx1024G/' ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-run.vmoptions
sed -i'' 's/-Xmx64G/-Xmx1024G/' ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-build.vmoptions

# install haystack 

python -m pip install git+https://github.com/antonisdim/Haystack

# create mutlifasta files for kraken2 and MALT and move them to the right folders

cat `cat 10_species.txt | paste -s -d ' '` > db_input_10_species.fasta.gz

cat `cat 100_species.txt | paste -s -d ' '` > db_input_100_species.fasta.gz 

cat `cat 500_species.txt | paste -s -d ' '` > db_input_500_species.fasta.gz 

cat `cat 1000_species.txt | paste -s -d ' '` > db_input_1000_species.fasta.gz

for i in ../rip_genome_cache/*; do cat $i/*fasta.gz >> db_input_5638_species.fasta.gz; done 

gunzip -d db_input_1000_species.fasta.gz 
gunzip -d db_input_100_species.fasta.gz
gunzip -d db_input_10_species.fasta.gz
gunzip -d db_input_500_species.fasta.gz
gunzip -d db_input_5638_species.fasta.gz

mkdir db_mutlifasta_inputs

mv *fasta db_mutlifasta_inputs/

# download mapping for MALT

mkdir mapping_files 
cd mapping_files
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Jul2020-2.db.zip 
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jul2020.db.zip
gunzip -c megan-map-Jul2020-2.db.zip > megan-map-Jul2020-2.db
cd ..

# make the database folders fpr Sigma

mkdir sigma_db_10_species
cd sigma_db_10_species/
for i in `less ../10_species.txt | cut -f 1-3 -d '/' | paste -s -d ' '`; do ln -s ../$i; done
cd ..

mkdir sigma_db_100_species
cd sigma_db_100_species/
for i in `less ../100_species.txt | cut -f 1-3 -d '/' | paste -s -d ' '`; do ln -s ../$i; done
cd ..

mkdir sigma_db_1000_species
cd sigma_db_1000_species/
for i in `less ../1000_species.txt | cut -f 1-3 -d '/' | paste -s -d ' '`; do ln -s ../$i; done
cd ..

mkdir sigma_db_500_species
cd sigma_db_500_species/
for i in `less ../500_species.txt | cut -f 1-3 -d '/' | paste -s -d ' '`; do ln -s ../$i; done
cd ..

mkdir sigma_db_5638_species
cd sigma_db_5638_species/
for i in `less ../5638_species.txt | cut -f 1-3 -d '/' | paste -s -d ' '`; do ln -s ../$i; done
cd ..

# run the performance tests

bash run_performance_test.sh &> performance_test.log 

