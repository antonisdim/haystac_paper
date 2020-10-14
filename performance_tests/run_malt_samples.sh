#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail


MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# run malt sample analysis performance test for samples of 1M reads against a db of 5638 species
/usr/bin/time -v malt-run \
	-t $MAX_CPU \
	-i ./inputs/input_1M_reads.fastq \
	-d index_new_5638_species \
	-m BlastN \
	-o .; \
	hops --mode me_po \
	--input input_1M_reads.rma6 \
	--output hops_input_1M_reads \
	--configFile configfile.txt

# run malt sample analysis performance test for samples of 10M reads against a db of 5638 species
/usr/bin/time -v malt-run \
	-t $MAX_CPU \
	-i ./inputs/input_10M_reads.fastq \
	-d index_new_5638_species \
	-m BlastN \
	-o .; \
	hops --mode me_po \
	--input input_10M_reads.rma6 \
	--output hops_input_10M_reads \
	--configFile configfile.txt

# run malt sample analysis performance test for samples of 100M reads against a db of 5638 species
/usr/bin/time -v malt-run \
	-t $MAX_CPU \
	-i ./inputs/input_100M_reads.fastq \
	-d index_new_5638_species \
	-m BlastN \
	-o .; \
	hops --mode me_po \
	--input input_100M_reads.rma6 \
	--output hops_input_100M_reads \
	--configFile configfile.txt

# run malt sample analysis performance test for samples of 10K reads against a db of 5638 species
/usr/bin/time -v malt-run \
	-t $MAX_CPU \
	-i ./inputs/input_10K_reads.fastq \
	-d index_new_5638_species \
	-m BlastN \
	-o .; \
	hops --mode me_po \
	--input input_10K_reads.rma6 \
	--output hops_input_10K_reads \
	--configFile configfile.txt

# run malt sample analysis performance test for samples of 100K reads against a db of 5638 species
/usr/bin/time -v malt-run \
	-t $MAX_CPU \
	-i ./inputs/input_100K_reads.fastq \
	-d index_new_5638_species \
	-m BlastN \
	-o .; \
	hops --mode me_po \
	--input input_100K_reads.rma6 \
	--output hops_input_100K_reads \
	--configFile configfile.txt

