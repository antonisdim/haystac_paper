time -v kraken2 --db db_kraken_5638_species --gzip-compressed --use-names --report input_1M.report --unclassified-out input_1M_unclassified.out --classified-out input_1M_classified.out --output input_1M.out ./inputs/input_1M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_5638_species -i input_1M.report -o input_1M.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_5638_species --gzip-compressed --use-names --report input_10K.report --unclassified-out input_10K_unclassified.out --classified-out input_10K_classified.out --output input_10K.out ./inputs/input_10K_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_5638_species -i input_10K.report -o input_10K.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_5638_species --gzip-compressed --use-names --report input_100K.report --unclassified-out input_100K_unclassified.out --classified-out input_100K_classified.out --output input_100K.out ./inputs/input_100K_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_5638_species -i input_100K.report -o input_100K.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_5638_species --gzip-compressed --use-names --report input_10M.report --unclassified-out input_10M_unclassified.out --classified-out input_10M_classified.out --output input_10M.out ./inputs/input_10M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_5638_species -i input_10M.report -o input_10M.bracken --threads `grep -c ^processor /proc/cpuinfo`

time -v kraken2 --db db_kraken_5638_species --gzip-compressed --use-names --report input_100M.report --unclassified-out input_100M_unclassified.out --classified-out input_100M_classified.out --output input_100M.out ./inputs/input_100M_reads.fastq.gz; Bracken-2.5/bracken -d db_kraken_5638_species -i input_100M.report -o input_100M.bracken --threads `grep -c ^processor /proc/cpuinfo`

