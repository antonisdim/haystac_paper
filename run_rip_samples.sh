conda activate rip

haystack config --bowtie2_threads 1

time -v haystack sample --sample-prefix input_10K_reads --fastq ./inputs/input_10K_reads.fastq.gz --output ./samples/input_10K_reads

time -v haystack analyse --mode abundances --database ./refseq_resp --sample ./samples/input_10K_reads --output ./refseq_analysis_output 

time -v haystack sample --sample-prefix input_100K_reads --fastq ./inputs/input_100K_reads.fastq.gz --output ./samples/input_100K_reads

time -v haystack analyse --mode abundances --database ./refseq_resp --sample ./samples/input_100K_reads --output ./refseq_analysis_output 

time -v haystack sample --sample-prefix input_100M_reads --fastq ./inputs/input_100M_reads.fastq.gz --output ./samples/input_100M_reads

time -v haystack analyse --mode abundances --database ./refseq_resp --sample ./samples/input_100M_reads --output ./refseq_analysis_output 

time -v haystack sample --sample-prefix input_1M_reads --fastq ./inputs/input_1M_reads.fastq.gz --output ./samples/input_1M_reads

time -v haystack analyse --mode abundances --database ./refseq_resp --sample ./samples/input_1M_reads --output ./refseq_analysis_output 

time -v haystack sample --sample-prefix input_10M_reads --fastq ./inputs/input_10M_reads.fastq.gz --output ./samples/input_10M_reads

time -v haystack analyse --mode abundances --database ./refseq_resp --sample ./samples/input_10M_reads --output ./refseq_analysis_output 

