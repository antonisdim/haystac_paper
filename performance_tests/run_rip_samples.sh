#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail


# run haystack sample analysis performance test for a sample of 10K reads against a db of 5638 species
/usr/bin/time -v haystack sample --sample-prefix input_10K_reads \
  --fastq ./inputs/input_10K_reads.fastq.gz --output ./samples/input_10K_reads

/usr/bin/time -v haystack analyse --mode abundances --database ./refseq_resp \
  --sample ./samples/input_10K_reads --output ./refseq_analysis_output

# run haystack sample analysis performance test for a sample of 100K reads against a db of 5638 species
/usr/bin/time -v haystack sample --sample-prefix input_100K_reads \
  --fastq ./inputs/input_100K_reads.fastq.gz --output ./samples/input_100K_reads

/usr/bin/time -v haystack analyse --mode abundances --database ./refseq_resp \
  --sample ./samples/input_100K_reads --output ./refseq_analysis_output

# run haystack sample analysis performance test for a sample of 100M reads against a db of 5638 species
/usr/bin/time -v haystack sample --sample-prefix input_100M_reads \
  --fastq ./inputs/input_100M_reads.fastq.gz --output ./samples/input_100M_reads

/usr/bin/time -v haystack analyse --mode abundances --database ./refseq_resp \
  --sample ./samples/input_100M_reads --output ./refseq_analysis_output

# run haystack sample analysis performance test for a sample of 1M reads against a db of 5638 species
/usr/bin/time -v haystack sample --sample-prefix input_1M_reads \
  --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

/usr/bin/time -v haystack analyse --mode abundances --database ./refseq_resp \
  --sample ./samples/input_1M_reads --output ./refseq_analysis_output

# run haystack sample analysis performance test for a sample of 10M reads against a db of 5638 species
/usr/bin/time -v haystack sample --sample-prefix input_10M_reads \
  --fastq ./inputs/input_10M_reads.fastq.gz --output ./samples/input_10M_reads

/usr/bin/time -v haystack analyse --mode abundances --database ./refseq_resp \
  --sample ./samples/input_10M_reads --output ./refseq_analysis_output

