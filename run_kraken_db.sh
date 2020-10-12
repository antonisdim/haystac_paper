time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_100_species.fasta --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
time -v kraken2-build --build --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
cp -r db_kraken db_kraken_100_species
rm -r db_kraken/library
rm -r db_kraken/database*
rm -r db_kraken/*k2d

time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_10_species.fasta --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
time -v kraken2-build --build --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
cp -r db_kraken db_kraken_10_species
rm -r db_kraken/library
rm -r db_kraken/database*
rm -r db_kraken/*k2d

time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_500_species.fasta --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
time -v kraken2-build --build --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
cp -r db_kraken db_kraken_500_species
rm -r db_kraken/library
rm -r db_kraken/database*
rm -r db_kraken/*k2d

time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_5638_species.fasta --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
time -v kraken2-build --build --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
cp -r db_kraken db_kraken_5638_species
rm -r db_kraken/library
rm -r db_kraken/database*
rm -r db_kraken/*k2d

time -v kraken2-build --add-to-library db_mutlifasta_inputs/db_input_1000_species.fasta --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
time -v kraken2-build --build --db db_kraken --threads `grep -c ^processor /proc/cpuinfo`
mv db_kraken db_kraken_1000_species

time -v Bracken-2.5/bracken-build -d db_kraken_100_species/
time -v Bracken-2.5/bracken-build -d db_kraken_10_species/ 
time -v Bracken-2.5/bracken-build -d db_kraken_1000_species/
time -v Bracken-2.5/bracken-build -d db_kraken_500_species/
time -v Bracken-2.5/bracken-build -d db_kraken_5638_species/
