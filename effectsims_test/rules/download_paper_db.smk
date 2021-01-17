#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


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
        "( haystac database --mode build --accessions entrez-selected-seqs.txt --refseq-rep --force-accessions "
        "--exclude-accessions AGIY02 --output paper_db --mem 64000 ) 2> {log}"
