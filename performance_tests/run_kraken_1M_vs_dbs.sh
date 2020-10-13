#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

time -v kraken2 --db db_kraken_100_species --gzip-compressed --use-names --report input_1M_100sp.report --unclassified-out input_1M_100sp_unclassified.out --classified-out input_1M_100sp_classified.out --output input_1M_100sp.out ./inputs/input_1M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_100_species -i input_1M_100sp.report -o input_1M_100sp.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_10_species --gzip-compressed --use-names --report input_1M_10sp.report --unclassified-out input_1M_10sp_unclassified.out --classified-out input_1M_10sp_classified.out --output input_1M_10sp.out ./inputs/input_1M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_10_species -i input_1M_10sp.report -o input_1M_10sp.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_1000_species --gzip-compressed --use-names --report input_1M_1000sp.report --unclassified-out input_1M_1000sp_unclassified.out --classified-out input_1M_1000sp_classified.out --output input_1M_1000sp.out ./inputs/input_1M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_1000_species -i input_1M_1000sp.report -o input_1M_1000sp.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_500_species --gzip-compressed --use-names --report input_1M_500sp.report --unclassified-out input_1M_500sp_unclassified.out --classified-out input_1M_500sp_classified.out --output input_1M_500sp.out ./inputs/input_1M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_500_species -i input_1M_500sp.report -o input_1M_500sp.bracken --threads `grep -c ^processor /proc/cpuinfo`
