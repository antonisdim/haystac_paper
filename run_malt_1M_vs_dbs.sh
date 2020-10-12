time -v malt-run -i ./inputs/input_1M_reads.fastq -d index_new_100_species -m BlastN -o input_1M_reads_100sp -t `grep -c ^processor /proc/cpuinfo`; hops --mode me_po --input input_1M_reads_100sp.rma6 --output hops_input_1M_reads_100sp --configFile configfile.txt

time -v malt-run -i ./inputs/input_1M_reads.fastq -d index_new_10_species -m BlastN -o input_1M_reads_10sp -t `grep -c ^processor /proc/cpuinfo`; hops --mode me_po --input input_1M_reads_10sp.rma6 --output hops_input_1M_reads_10sp --configFile configfile.txt

time -v malt-run -i ./inputs/input_1M_reads.fastq -d index_new_500_species -m BlastN -o input_1M_reads_500sp -t `grep -c ^processor /proc/cpuinfo`; hops --mode me_po --input input_1M_reads_500sp.rma6 --output hops_input_1M_reads_500sp --configFile configfile.txt

time -v malt-run -i ./inputs/input_1M_reads.fastq -d index_new_1000_species -m BlastN -o input_1M_reads_1000sp -t `grep -c ^processor /proc/cpuinfo`; hops --mode me_po --input input_1M_reads_1000sp.rma6 --output hops_input_1M_reads_1000sp --configFile configfile.txt

