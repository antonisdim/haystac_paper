time -v malt-build -t `grep -c ^processor /proc/cpuinfo` -s DNA -i db_mutlifasta_inputs/db_input_100_species.fasta -d index_new/ -a2taxonomy mapping_files/megan-map-Oct2019.db
mv index_new index_new_100_species
time -v malt-build -t `grep -c ^processor /proc/cpuinfo` -s DNA -i db_mutlifasta_inputs/db_input_10_species.fasta -d index_new/ -a2taxonomy mapping_files/megan-map-Oct2019.db
mv index_new index_new_10_species
time -v malt-build -t `grep -c ^processor /proc/cpuinfo` -s DNA -i db_mutlifasta_inputs/db_input_500_species.fasta -d index_new/ -a2taxonomy mapping_files/megan-map-Oct2019.db
mv index_new index_new_500_species
time -v malt-build -t `grep -c ^processor /proc/cpuinfo` -s DNA -i db_mutlifasta_inputs/db_input_5638_species.fasta -d index_new/ -a2taxonomy mapping_files/megan-map-Oct2019.db
mv index_new mv index_new index_new_5638_species
time -v malt-build -t `grep -c ^processor /proc/cpuinfo` -s DNA -i db_mutlifasta_inputs/db_input_1000_species.fasta -d index_new/ -a2taxonomy mapping_files/megan-map-Oct2019.db
mv index_new mv index_new index_new_1000_species
