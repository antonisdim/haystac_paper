#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

# bash strict mode
set -euo pipefail

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

# run kraken building db performance test for 100 species
/usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_100_species
/usr/bin/time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_100_species.fasta \
  --db db_kraken_100_species --threads $MAX_CPU
/usr/bin/time -v kraken2-build --build --db db_kraken_100_species --threads $MAX_CPU

# run kraken building db performance test for 10 species
/usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_10_species
/usr/bin/time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_10_species.fasta \
  --db db_kraken_10_species --threads $MAX_CPU
/usr/bin/time -v kraken2-build --build --db db_kraken_10_species --threads $MAX_CPU

# run kraken building db performance test for 500 species
/usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_500_species
/usr/bin/time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_500_species.fasta \
  --db db_kraken_500_species --threads $MAX_CPU
/usr/bin/time -v kraken2-build --build --db db_kraken_500_species --threads $MAX_CPU

# run kraken building db performance test for 5638 species
/usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_5638_species
/usr/bin/time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_5638_species.fasta \
  --db db_kraken_5638_species --threads $MAX_CPU
/usr/bin/time -v kraken2-build --build --db db_kraken_5638_species --threads $MAX_CPU

# run kraken building db performance test for 1000 species
/usr/bin/time -v kraken2-build --download-taxonomy --db db_kraken_1000_species
/usr/bin/time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_1000_species.fasta \
  --db db_kraken_1000_species --threads $MAX_CPU
/usr/bin/time -v kraken2-build --build --db db_kraken_1000_species --threads $MAX_CPU

# run bracken building db performance test for 100 species
/usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_100_species/
# run bracken building db performance test for 10 species
/usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_10_species/
# run bracken building db performance test for 1000 species
/usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_1000_species/
# run bracken building db performance test for 500 species
/usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_500_species/
# run bracken building db performance test for 5638 species
/usr/bin/time -v Bracken-2.5/bracken-build -d db_kraken_5638_species/
