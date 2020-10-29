#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# download a large HMP sample
mkdir -p inputs
wget --quiet http://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v2/SRS078671.tar.bz2 \
  -O inputs/SRS078671.tar.bz2

# decompress the downloaded file
tar -xvjf inputs/SRS078671.tar.bz2

# concatenate the files into one big file
cat inputs/SRS078671/SRS078671.denovo_duplicates_marked.trimmed.1.fastq \
  inputs/SRS078671/SRS078671.denovo_duplicates_marked.trimmed.singleton.fastq > \
  inputs/SRS078671_aggregate.fastq

# subsample the file into smaller samples
read_counts=(input_10K input_100K input_1M input_10M input_100M)
size=1000

for read_count in "${read_counts[@]}"; do
  echo "Subsampling for ${read_count}"
  size=$((size*10))
  seqtk sample inputs/SRS078671_aggregate.fastq $size > inputs/"${read_count}"_reads.fastq
done

# compress the files
for read_count in "${read_counts[@]}"; do
  echo "Compressing the subsampled file for ${read_count}"
  gzip -c inputs/"${read_count}"_reads.fastq > inputs/"${read_count}"_reads.fastq.gz
done

# removing files that are not going to be used
rm inputs/SRS078671_aggregate.fastq
rm -r inputs/SRS078671
