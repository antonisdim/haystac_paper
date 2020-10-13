#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# run malt sample analysis performance test against dbs of various sizes

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

time -v malt-run \
	-i ./inputs/input_1M_reads.fastq \
	-d index_new_100_species \
	-m BlastN \
	-o input_1M_reads_100sp \
	-t $MAX_CPU; \
	hops --mode me_po \
	--input input_1M_reads_100sp.rma6 \
	--output hops_input_1M_reads_100sp \
	--configFile configfile.txt

time -v malt-run \
	-i ./inputs/input_1M_reads.fastq \
	-d index_new_10_species \
	-m BlastN \
	-o input_1M_reads_10sp \
	-t $MAX_CPU; \
	hops --mode me_po \
	--input input_1M_reads_10sp.rma6 \
	--output hops_input_1M_reads_10sp \
	--configFile configfile.txt

time -v malt-run \
	-i ./inputs/input_1M_reads.fastq \
	-d index_new_500_species \
	-m BlastN \
	-o input_1M_reads_500sp \	
	-t $MAX_CPU; \
	hops --mode me_po \
	--input input_1M_reads_500sp.rma6 \
	--output hops_input_1M_reads_500sp \
	--configFile configfile.txt

time -v malt-run \
	-i ./inputs/input_1M_reads.fastq \
	-d index_new_1000_species \
	-m BlastN \
	-o input_1M_reads_1000sp \
	-t $MAX_CPU; \
	hops --mode me_po \
	--input input_1M_reads_1000sp.rma6 \
	--output hops_input_1M_reads_1000sp \
	--configFile configfile.txt

