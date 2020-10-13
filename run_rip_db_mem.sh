# run test for db building with haystack with memory limit  

MEM_LIMIT_TEST = 8000

haystack config --email antonisdim41@gmail.com
haystack config --genome-cache-folder ../rip_genome_cache/
haystack config --batchsize 5
haystack config --mismatch-probability 0.05
haystack config --bowtie2-threads 5
haystack config --bowtie2-scaling 15

rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_10_species_input.txt -o ./rip_db_10_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_100_species_input.txt -o ./rip_db_100_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_500_species_input.txt -o ./rip_db_500_species_input --mem $MEM_LIMIT_TEST
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database --query '("Yersinia"[Organism] OR "Haemophilus"[Organism] OR "Klebsiella"[Organism] OR "Bordetella"[Organism] OR "Streptococcus"[Organism]) AND "complete genome"[All Fields] AND refseq[filter]' \
	--output ./refseq_resp --mem $MEM_LIMIT_TEST \
	--refseq-rep True
rm ../rip_genome_cache/*/*.bt2l
time -v haystack database -a rip_db_1000_species_input.txt -o ./rip_db_1000_species_input --mem $MEM_LIMIT_TEST
