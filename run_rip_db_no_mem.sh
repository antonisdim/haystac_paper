#run 
conda activate rip 

haystack config --email antonis1909@hotmail.com
haystack config --genome-cache-folder `readlink -f ../rip_genome_cache/`
haystack config --bowtie2_scaling 15
haystack config --bowtie2_threads 5

for i in ../rip_genome_cache/*; do rm $i/*.bt2l; done
time -v haystack database -a rip_db_10_species_input.txt -o ./rip_db_10_species_input_no_mem
for i in ../rip_genome_cache/*; do rm $i/*.bt2l; done
time -v haystack database -a rip_db_100_species_input.txt -o ./rip_db_100_species_input_no_mem
for i in ../rip_genome_cache/*; do rm $i/*.bt2l; done
time -v haystack database -a rip_db_500_species_input.txt -o ./rip_db_500_species_input_no_mem
for i in ../rip_genome_cache/*; do rm $i/*.bt2l; done
time -v haystack database --query '("Yersinia"[Organism] OR "Haemophilus"[Organism] OR "Klebsiella"[Organism] OR "Bordetella"[Organism] OR "Streptococcus"[Organism]) AND "complete genome"[All Fields] AND refseq[filter]' --output ./refseq_resp_no_mem --refseq-rep True
for i in ../rip_genome_cache/*; do rm $i/*.bt2l; done
time -v haystack database -a rip_db_1000_species_input.txt -o ./rip_db_1000_species_input_no_mem
