#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import yaml
from os.path import expanduser
import pandas as pd


HAYSTAC_CONFIG = expanduser("~/.haystac/config.yaml")

with open(HAYSTAC_CONFIG, "r") as fin:
    conf = yaml.safe_load(fin)

CACHE = conf['cache']


rule prepare_garg_dirs:
    input:
        "paper_db/db_taxa_accessions.tsv",
    output:
        "garg_sims/{taxon}/endo/{accession}.fasta"
    message:
        "Preparing the gargammel inputs for species {wildcards.taxon} with accession {wildcards.accession}."
    shell:
        "mkdir -p garg_sims/{wildcards.taxon}/endo garg_sims/{wildcards.taxon}/bact garg_sims/{wildcards.taxon}/cont;"
        "gzip -c -d {CACHE}/ncbi/{wildcards.taxon}/{wildcards.accession}.fasta.gz > {output}"


def get_tax_acc_pairs(_):
    """Function to get the pairs of tax and acc for the gargammel sims"""

    db = checkpoints.download_paper_db.get()

    tax_acc = pd.read_csv(db.output.db_list, sep = '\t', names=['Taxon', 'Accession'])
    genbank_plasmids = pd.read_csv(db.output.gen_plasmids, sep = '\t')
    refseq_plasmids = pd.read_csv(db.output.ref_plasmids, sep = '\t')

    # drop plasmids - gargammel can only do one fasta per folder

    cond_gen = tax_acc['Accession'].isin(genbank_plasmids['AccessionVersion'])
    tax_acc.drop(tax_acc[cond_gen].index, inplace = True)

    cond_ref = tax_acc['Accession'].isin(refseq_plasmids['AccessionVersion'])
    tax_acc.drop(tax_acc[cond_ref].index, inplace = True)

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row['Taxon']
        accession = row['Accession']
        paths.append(f"garg_sims/{taxon}/endo/{accession}.fasta")

    return paths


rule all_dirs_prepared:
    input:
        get_tax_acc_pairs
    output:
        "all_garg_inputs_prepared.done"
    message:
        "All gargammel inputs have been prepared."
    shell:
        "touch {output}"


