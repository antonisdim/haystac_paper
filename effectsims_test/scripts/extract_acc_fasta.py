#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import sys
from Bio import SeqIO


def extract_acc_fasta(fasta_files, output_file):
    """Function to extract all the accessions in the fasta files of the simulated genomes"""

    print(
        "Creating the tsv file with all the simulated taxa and all their accessions ...",
        file=sys.stderr,
    )

    with open(output_file, "w") as fout:

        for fasta_file in fasta_files:

            print("Adding accessions from ", fasta_file, file=sys.stderr)

            # path format garg_sims/Piscicoccus_intestinalis/endo/BCNQ01.fasta
            taxon = fasta_file.split("/")[1]
            records = SeqIO.parse(fasta_file, "fasta")
            rec_list = [rec.id for rec in records]

            for acc in rec_list:
                print(taxon, acc, sep="\t", file=fout)

    print("TSV file ready", file=sys.stderr)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    extract_acc_fasta(fasta_files=snakemake.input, output_file=snakemake.output[0])
