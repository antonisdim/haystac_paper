#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# create a new conda environment to run all the test_sets
eval "$(conda shell.bash hook)"
conda env remove --name small_test
conda env create -f environment.yaml -n small_test
conda activate small_test

# bash strict mode (set after conda is activated, as the start-up scripts are not strict safe)
set -euo pipefail

# force bash to use the `time` command and not it's own native implementation
time=$(which time)

# set some config options for haystack
haystack config \
  --api-key 4a7ceec35df93baf39e6f747e08f69f39909 \
  --cache ./haystack_genome_cache_small/

# fetch all the genomes in the representative RefSeq database
$time -v haystack database \
  --mode fetch \
  --accessions ./haystack_configs/haystack_db_100_species_refseq_input.txt \
  --output ./genome_fetch_small/ \
  --force-accessions

# download a large HMP sample
mkdir -p inputs
wget --quiet http://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v2/SRS078671.tar.bz2 \
  -O inputs/SRS078671.tar.bz2

# decompress the downloaded file
tar -xvjf inputs/SRS078671.tar.bz2 -C inputs/

# concatenate the files into one big file
cat inputs/SRS078671/SRS078671.denovo_duplicates_marked.trimmed.1.fastq \
  inputs/SRS078671/SRS078671.denovo_duplicates_marked.trimmed.singleton.fastq > \
  inputs/SRS078671_aggregate.fastq

# subsample the file into smaller samples
read_counts=(input_10K input_100K)
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

echo "DONE!"
