#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

REAL_HUMAN = (
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR218/007/ERR2181127/ERR2181127.fastq.gz"
)


checkpoint download_paper_db:
    input:
        genome_list="entrez-selected-seqs.txt",
    log:
        "download_paper_db.log",
    output:
        db_list="paper_db/db_taxa_accessions.tsv",
        gen_plasmids="paper_db/entrez/genbank-plasmids.tsv",
        ref_plasmids="paper_db/entrez/refseq-plasmids.tsv",
    message:
        "Building the prok representative RefSeq with an additional NCBI query."
    shell:
        "( haystac database --mode build --accessions entrez-selected-seqs.txt --refseq-rep prokaryote_rep "
        "--force-accessions --exclude-accessions AGIY02 --output paper_db --mem 64000 ) 2> {log}"


checkpoint download_sims_accs:
    input:
        genome_list="garg_sim_accessions.tsv",
    log:
        "sim_acc_paper_db.log",
    output:
        sim_list="sim_accs/db_taxa_accessions.tsv",
    message:
        "Downloading the gargammel simluation accessions."
    shell:
        "( haystac database --mode fetch --accessions {input} --force-accessions "
        "--exclude-accessions AGIY02 --output sim_accs --mem 64000 ) 2> {log}"


checkpoint download_spike_accs:
    input:
        genome_list="garg_spike_species.tsv",
    log:
        "spike_acc_paper_db.log",
    output:
        sim_list="spike_accs/db_taxa_accessions.tsv",
    message:
        "Downloading the gargammel simluation spike accessions."
    shell:
        "( haystac database --mode fetch --accessions {input} --force-accessions "
        "--exclude-accessions AGIY02 --output spike_accs --mem 64000 ) 2> {log}"


checkpoint download_real_human:
    log:
        "real_human/ERR2181127.log",
    output:
        "real_human/ERR2181127.fastq.gz",
    message:
        "Downloading a real human genome from Lipson et al 2017."
    shell:
        "wget -O {output} {REAL_HUMAN} 2> {log}"
