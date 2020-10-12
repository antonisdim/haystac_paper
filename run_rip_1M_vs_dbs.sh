conda activate rip

haystack config --bowtie2_threads 1

rm -r ./samples/input_1M_reads

time -v haystack sample --sample-prefix input_1M_reads --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

time -v haystack analyse --mode abundances --database ./rip_db_10_species_input --sample ./samples/input_1M_reads --output ./rip_db_10_species_analysis_output 

rm -r ./samples/input_1M_reads

time -v haystack sample --sample-prefix input_1M_reads --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

time -v haystack analyse --mode abundances --database ./rip_db_100_species_input --sample ./samples/input_1M_reads --output ./rip_db_100_species_analysis_output 

rm -r ./samples/input_1M_reads

time -v haystack sample --sample-prefix input_1M_reads --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

time -v haystack analyse --mode abundances --database ./rip_db_500_species_input --sample ./samples/input_1M_reads --output ./rip_db_500_species_analysis_output 

rm -r ./samples/input_1M_reads

time -v haystack sample --sample-prefix input_1M_reads --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

time -v haystack analyse --mode abundances --database ./rip_db_1000_species_input --sample ./samples/input_1M_reads --output ./rip_db_1000_species_analysis_output 
