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

CACHE = conf["cache"]


rule prepare_garg_dirs:
    input:
        "paper_db/db_taxa_accessions.tsv",
    output:
        "garg_sims/{taxon}/endo/{accession}.fasta",
    message:
        "Preparing the gargammel inputs for species {wildcards.taxon} with accession {wildcards.accession}."
    shell:
        "mkdir -p garg_sims/{wildcards.taxon}/endo garg_sims/{wildcards.taxon}/bact garg_sims/{wildcards.taxon}/cont;"
        "gzip -c -d {CACHE}/ncbi/{wildcards.taxon}/{wildcards.accession}.fasta.gz > {output}"


def get_tax_acc_pairs(_):
    """Function to get the pairs of tax and acc for the gargammel sims"""

    db = checkpoints.download_paper_db.get()

    tax_acc = pd.read_csv(db.output.db_list, sep="\t", names=["Taxon", "Accession"])
    genbank_plasmids = pd.read_csv(db.output.gen_plasmids, sep="\t")
    refseq_plasmids = pd.read_csv(db.output.ref_plasmids, sep="\t")

    # drop plasmids - gargammel can only do one fasta per folder

    cond_gen = tax_acc["Accession"].isin(genbank_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_gen].index, inplace=True)

    cond_ref = tax_acc["Accession"].isin(refseq_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_ref].index, inplace=True)

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row["Taxon"]
        accession = row["Accession"]
        paths.append(f"garg_sims/{taxon}/endo/{accession}.fasta")

    return paths


rule all_dirs_prepared:
    input:
        get_tax_acc_pairs,
    output:
        "all_garg_inputs_prepared.done",
    message:
        "All gargammel inputs have been prepared."
    shell:
        "touch {output}"


rule simulate_gargammel:
    input:
        "all_garg_inputs_prepared.done",
    log:
        "garg_sims/{taxon}.log",
    output:
        "garg_sims/{taxon}_s.fq.gz",
    message:
        "Simulating 1000 reads from taxon {wildcards.taxon}."
    conda:
        "../envs/gargammel.yaml"
    shell:
        "( gargammel --comp 0,0,1 -n 1000 -l 100 single -rl 100 -se -ss HS25 garg_sims/{wildcards.taxon}/ "
        "-o garg_sims/{wildcards.taxon} ) &> {log}"


def get_tax_frags(_):
    """Function to get the sim frags per taxon"""
    # todo I know it is duplicated, but is it worth it putting it in a utilities file it's only one function ?

    db = checkpoints.download_paper_db.get()

    tax_acc = pd.read_csv(db.output.db_list, sep="\t", names=["Taxon", "Accession"])
    genbank_plasmids = pd.read_csv(db.output.gen_plasmids, sep="\t")
    refseq_plasmids = pd.read_csv(db.output.ref_plasmids, sep="\t")

    # drop plasmids - gargammel can only do one fasta per folder

    cond_gen = tax_acc["Accession"].isin(genbank_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_gen].index, inplace=True)

    cond_ref = tax_acc["Accession"].isin(refseq_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_ref].index, inplace=True)

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row["Taxon"]
        paths.append(f"garg_sims/{taxon}_s.fq.gz")

    return paths


rule create_meta_sim:
    input:
        get_tax_frags,
    log:
        "raw_samples/effectsim_lib.log",
    output:
        "raw_samples/effectsim_lib.fastq.gz",
    message:
        "Creating a library of approx 5600 species, with 100 bp read len and 0 deamination rate."
    script:
        "../scripts/copy.py"


rule vary_length_damage:
    input:
        "raw_samples/effectsim_lib.fastq.gz",
    log:
        "raw_samples/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        "raw_samples/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.fastq.gz",
    conda:
        "../envs/gargammel.yaml"
    message:
        "Changing the read length to {wildcards.readlen} bp and the deamination rate to {wildcards.damage}."
    params:
        bptotrim=lambda wildcards: 100 - int(wildcards.readlen),
        # art_len=lambda wildcards: 120 if int(wildcards.readlen) != 120 else 108,
    shell:
        "(seqtk trimfq -e {params.bptotrim} {input} > "
        "raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq;"
        "seqtk seq -a raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq > "
        "raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta;"
        "deamSim -damage 0.03,{wildcards.overhang},0.01,{wildcards.damage} "
        "raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta > "
        "raw_samples/deamsim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta;"
        "art_illumina -ss HS25 -amp -na "
        "-i raw_samples/deamsim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta "
        "-l {wildcards.readlen} -c 1 -qs 0 -qs2 0 "
        "-o raw_samples/illumina_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l; "
        "bgzip --force "
        "--stdout raw_samples/illumina_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fq > {output}; "
        "rm raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq; "
        "rm raw_samples/seq_trim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta; "
        "rm raw_samples/deamsim_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fasta; "
        "rm raw_samples/illumina_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fq) 2> {log}"


# simulation-{id}_readlen-{readlen}_damage-{damage}.fq.gz
