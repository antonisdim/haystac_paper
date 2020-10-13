#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run malt building db performance test for samples of various sizes

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

time -v malt-build \
	-t $MAX_CPU \
	-s DNA \
	-i db_mutlifasta_inputs/db_input_100_species.fasta \
	-d index_new_100_species/ \
	-a2taxonomy mapping_files/megan-map-Jul2020-2.db 

time -v malt-build \
	-t $MAX_CPU \
	-s DNA \
	-i db_mutlifasta_inputs/db_input_10_species.fasta \
	-d index_new_10_species/ \
	-a2taxonomy mapping_files/megan-map-Jul2020-2.db 

time -v malt-build \
	-t $MAX_CPU \
	-s DNA \
	-i db_mutlifasta_inputs/db_input_500_species.fasta \
	-d index_new_500_species/ \
	-a2taxonomy mapping_files/megan-map-Jul2020-2.db

time -v malt-build \
	-t $MAX_CPU \
	-s DNA \
	-i db_mutlifasta_inputs/db_input_5638_species.fasta \
	-d index_new_5638_species/ \
	-a2taxonomy mapping_files/megan-map-Jul2020-2.db

time -v malt-build \
	-t $MAX_CPU \
	-s DNA \
	-i db_mutlifasta_inputs/db_input_1000_species.fasta \
	-d index_new_1000_species/ \
	-a2taxonomy mapping_files/megan-map-Jul2020-2.db
