#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# run haystack building db performance test with no mem limit

rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_10_species_input.txt -o ./rip_db_10_species_input_no_mem
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_100_species_input.txt -o ./rip_db_100_species_input_no_mem
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_500_species_input.txt -o ./rip_db_500_species_input_no_mem
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database --query '("Yersinia"[Organism] OR "Haemophilus"[Organism] OR "Klebsiella"[Organism] OR "Bordetella"[Organism] OR "Streptococcus"[Organism]) AND "complete genome"[All Fields] AND refseq[filter]' --output ./refseq_resp_no_mem --refseq-rep True
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_1000_species_input.txt -o ./rip_db_1000_species_input_no_mem
