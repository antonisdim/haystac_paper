#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
import sys
from Bio import SeqIO


def extract_acc_fasta(fastq_files, output_file):
    """Function to extract all the accessions in the fasta files of the simulated genomes"""

    print(
        "Creating the tsv file with all the simulated taxa and all their accessions ...",
        file=sys.stderr,
    )

    with open(output_file, "w") as fout:

        for fastq_file in fastq_files:

            print("Adding read from ", fastq_file, file=sys.stderr)

            # path format garg_sims/Corynebacterium_diphtheriae_s.fq.gz
            taxon = fastq_file.split("/")[1].replace("_s.fq.gz", "")
            with gzip.open(fastq_file, "rt") as fastq_handle:
                records = SeqIO.parse(fastq_handle, "fastq")
                rec_list = [rec.id for rec in records]

            for acc in rec_list:
                print(taxon, acc, sep="\t", file=fout)

    print("TSV file ready", file=sys.stderr)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    extract_acc_fasta(fastq_files=snakemake.input, output_file=snakemake.output[0])
