#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import shutil
import sys
import gzip
from Bio import bgzf


def bowtie2_multifasta(fasta_files, output_file):
    print(
        "Creating the bowtie2 multifasta file for the filtering ...", file=sys.stderr,
    )

    # TODO this is super slow, see https://sites.google.com/site/tfsidc/linux-tricks/processing-a-large-number-of-files
    with bgzf.open(output_file, "wt") as fout:
        for fasta_file in fasta_files:
            print("Adding reads from ", fasta_file, file=sys.stderr)
            with gzip.open(fasta_file, "rb") as fin:
                shutil.copyfileobj(fin, fout)

    print("Multi fastq file created", file=sys.stderr)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    bowtie2_multifasta(fasta_files=snakemake.input, output_file=snakemake.output[0])
