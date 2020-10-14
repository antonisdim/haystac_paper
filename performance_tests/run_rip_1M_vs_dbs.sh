#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

# run haystack sample analysis performance test for a sample of 1M reads against a db of 10 species
rm -r ./samples/input_1M_reads

/usr/bin/time -v haystack sample --sample-prefix input_1M_reads \
  --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

/usr/bin/time -v haystack analyse \
	--mode abundances \
	--database ./rip_db_10_species_input \
	--sample ./samples/input_1M_reads \
	--output ./rip_db_10_species_analysis_output 

# run haystack sample analysis performance test for a sample of 1M reads against a db of 100 species
rm -r ./samples/input_1M_reads

/usr/bin/time -v haystack sample --sample-prefix input_1M_reads \
  --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

/usr/bin/time -v haystack analyse \
	--mode abundances \
	--database ./rip_db_100_species_input \
	--sample ./samples/input_1M_reads \
	--output ./rip_db_100_species_analysis_output 

# run haystack sample analysis performance test for a sample of 1M reads against a db of 500 species
rm -r ./samples/input_1M_reads

/usr/bin/time -v haystack sample --sample-prefix input_1M_reads \
  --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

/usr/bin/time -v haystack analyse \
	--mode abundances \
	--database ./rip_db_500_species_input \
	--sample ./samples/input_1M_reads \
	--output ./rip_db_500_species_analysis_output 

# run haystack sample analysis performance test for a sample of 1M reads against a db of 1000 species
rm -r ./samples/input_1M_reads

/usr/bin/time -v haystack sample --sample-prefix input_1M_reads \
  --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

/usr/bin/time -v haystack analyse \
	--mode abundances \
	--database ./rip_db_1000_species_input \
	--sample ./samples/input_1M_reads \
	--output ./rip_db_1000_species_analysis_output 
