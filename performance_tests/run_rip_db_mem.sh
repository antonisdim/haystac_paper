#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail



MEM_LIMIT_TEST=8000

haystack config \
  --email antonisdim41@gmail.com \
  --genome-cache-folder ../rip_genome_cache/

# run test for db building with haystack with memory limit for a db of 10 species
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a ./haystack_configs/rip_db_10_species_input.txt \
  -o ./rip_db_10_species_input --mem $MEM_LIMIT_TEST

# run test for db building with haystack with memory limit for a db of 100 species
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a ./haystack_configs/rip_db_100_species_input.txt \
  -o ./rip_db_100_species_input --mem $MEM_LIMIT_TEST

# run test for db building with haystack with memory limit for a db of 500 species
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a ./haystack_configs/rip_db_500_species_input.txt \
-o ./rip_db_500_species_input --mem $MEM_LIMIT_TEST

# run test for db building with haystack with memory limit for a db of 5638 species
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database --query '("Yersinia"[Organism] OR "Haemophilus"[Organism] OR "Klebsiella"[Organism] OR "Bordetella"[Organism] OR "Streptococcus"[Organism]) AND "complete genome"[All Fields] AND refseq[filter]' \
	--output ./refseq_resp --mem $MEM_LIMIT_TEST \
	--refseq-rep True

# run test for db building with haystack with memory limit for a db of 1000 species
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a ./haystack_configs/rip_db_1000_species_input.txt \
-o ./rip_db_1000_species_input --mem $MEM_LIMIT_TEST
