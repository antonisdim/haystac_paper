#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# run kraken sample analysis performance test for a sample with 1M reads
/usr/bin/time -v kraken2 \
	--db db_kraken_5638_species \
	--gzip-compressed \
	--use-names \
	--report input_1M.report \
	--unclassified-out input_1M_unclassified.out \
	--classified-out input_1M_classified.out \
	--threads $MAX_CPU \
	--output input_1M.out ./inputs/input_1M_reads.fastq.gz; \
	Bracken-2.5/bracken -d db_kraken_5638_species \
	-i input_1M.report \
	-o input_1M.bracken \
	--threads $MAX_CPU

# run kraken sample analysis performance test for a sample with 10K reads
/usr/bin/time -v kraken2 \
	--db db_kraken_5638_species \
	--gzip-compressed \
	--use-names \
	--report input_10K.report \
	--unclassified-out input_10K_unclassified.out \
	--classified-out input_10K_classified.out \
	--threads $MAX_CPU \
	--output input_10K.out ./inputs/input_10K_reads.fastq.gz; \
	Bracken-2.5/bracken -d db_kraken_5638_species \
	-i input_10K.report \
	-o input_10K.bracken \
	--threads $MAX_CPU

# run kraken sample analysis performance test for a sample with 100K reads
/usr/bin/time -v kraken2 \
	--db db_kraken_5638_species \
	--gzip-compressed \
	--use-names \
	--report input_100K.report \
	--unclassified-out input_100K_unclassified.out \
	--classified-out input_100K_classified.out \
	--threads $MAX_CPU \
	--output input_100K.out ./inputs/input_100K_reads.fastq.gz; \
	Bracken-2.5/bracken -d db_kraken_5638_species \
	-i input_100K.report \
	-o input_100K.bracken \
	--threads $MAX_CPU

# run kraken sample analysis performance test for a sample with 10M reads
/usr/bin/time -v kraken2 \
	--db db_kraken_5638_species \
	--gzip-compressed \
	--use-names \
	--report input_10M.report \
	--unclassified-out input_10M_unclassified.out \
	--classified-out input_10M_classified.out \
	--threads $MAX_CPU \
	--output input_10M.out ./inputs/input_10M_reads.fastq.gz; \
	Bracken-2.5/bracken -d db_kraken_5638_species \
	-i input_10M.report \
	-o input_10M.bracken \
	--threads $MAX_CPU

# run kraken sample analysis performance test for a sample with 100M reads
/usr/bin/time -v kraken2 \
	--db db_kraken_5638_species \
	--gzip-compressed \
	--use-names \
	--report input_100M.report \
	--unclassified-out input_100M_unclassified.out \
	--classified-out input_100M_classified.out \
	--threads $MAX_CPU \
	--output input_100M.out ./inputs/input_100M_reads.fastq.gz; \
	Bracken-2.5/bracken -d db_kraken_5638_species \
	-i input_100M.report \
	-o input_100M.bracken \
	--threads $MAX_CPU
