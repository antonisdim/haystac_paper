time -v malt-run -t `grep -c ^processor /proc/cpuinfo` -i ./inputs/input_1M_reads.fastq -d index_new_5638_species -m BlastN -o .; hops --mode me_po --input input_1M_reads.rma6 --output hops_input_1M_reads --configFile configfile.txt

time -v malt-run -t `grep -c ^processor /proc/cpuinfo` -i ./inputs/input_10M_reads.fastq -d index_new_5638_species -m BlastN -o .; hops --mode me_po --input input_10M_reads.rma6 --output hops_input_10M_reads --configFile configfile.txt

time -v malt-run -t `grep -c ^processor /proc/cpuinfo` -i ./inputs/input_100M_reads.fastq -d index_new_5638_species -m BlastN -o .; hops --mode me_po --input input_100M_reads.rma6 --output hops_input_100M_reads --configFile configfile.txt

time -v malt-run -t `grep -c ^processor /proc/cpuinfo` -i ./inputs/input_10K_reads.fastq -d index_new_5638_species -m BlastN -o .; hops --mode me_po --input input_10K_reads.rma6 --output hops_input_10K_reads --configFile configfile.txt

time -v malt-run -t `grep -c ^processor /proc/cpuinfo` -i ./inputs/input_100K_reads.fastq -d index_new_5638_species -m BlastN -o .; hops --mode me_po --input input_100K_reads.rma6 --output hops_input_100K_reads --configFile configfile.txt

