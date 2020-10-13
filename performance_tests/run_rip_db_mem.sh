#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run test for db building with haystack with memory limit

MEM_LIMIT_TEST=8000

haystack config \
  --email antonisdim41@gmail.com \
  --genome-cache-folder ../rip_genome_cache/

rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a rip_db_10_species_input.txt -o ./rip_db_10_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a rip_db_100_species_input.txt -o ./rip_db_100_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a rip_db_500_species_input.txt -o ./rip_db_500_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database --query '("Yersinia"[Organism] OR "Haemophilus"[Organism] OR "Klebsiella"[Organism] OR "Bordetella"[Organism] OR "Streptococcus"[Organism]) AND "complete genome"[All Fields] AND refseq[filter]' \
	--output ./refseq_resp --mem $MEM_LIMIT_TEST \
	--refseq-rep True
rm ../rip_genome_cache/*/*.bt2l
/usr/bin/time -v haystack database -a rip_db_1000_species_input.txt -o ./rip_db_1000_species_input --mem $MEM_LIMIT_TEST
