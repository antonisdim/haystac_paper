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
        "sim_accs/db_taxa_accessions.tsv",
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
    sims = checkpoints.download_sims_accs.get()

    tax_acc = pd.read_csv(sims.output.sim_list, sep="\t", names=["Taxon", "Accession"])
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
        "( gargammel --comp 0,0,1 -n 1000 -l 120 single -rl 120 -se -ss HS25 garg_sims/{wildcards.taxon}/ "
        "-o garg_sims/{wildcards.taxon} ) &> {log}"


def get_tax_frags(wildcards):
    """Function to get the sim frags per taxon"""
    # todo I know it is duplicated, but is it worth it putting it in a utilities file it's only one function ?

    db = checkpoints.download_paper_db.get()
    sims = checkpoints.download_sims_accs.get()

    tax_acc = pd.read_csv(sims.output.sim_list, sep="\t", names=["Taxon", "Accession"])
    genbank_plasmids = pd.read_csv(db.output.gen_plasmids, sep="\t")
    refseq_plasmids = pd.read_csv(db.output.ref_plasmids, sep="\t")

    # drop plasmids - gargammel can only do one fasta per folder

    cond_gen = tax_acc["Accession"].isin(genbank_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_gen].index, inplace=True)

    cond_ref = tax_acc["Accession"].isin(refseq_plasmids["AccessionVersion"])
    tax_acc.drop(tax_acc[cond_ref].index, inplace=True)

    if int(wildcards.species) != 5652:
        tax_acc = tax_acc.sample(int(wildcards.species))

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row["Taxon"]
        paths.append(f"garg_sims/{taxon}_s.fq.gz")

    return paths


rule create_meta_sim:
    input:
        get_tax_frags,
    log:
        "raw_samples/effectsim_lib_{species}sp.log",
    output:
        "raw_samples/effectsim_lib_{species}sp.fastq.gz",
    message:
        "Creating a library of {wildcards.species} species, with 120 bp read len and 0 deamination rate."
    script:
        "../scripts/copy.py"


rule seqtk:
    input:
        "raw_samples/effectsim_lib_{species}sp.fastq.gz",
    log:
        "raw_samples/seq_trim_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        temp("raw_samples/seq_trim_{species}sp_{readlen}bp_{damage}d_{overhang}l.fasta"),
    conda:
        "../envs/gargammel.yaml"
    message:
        "Changing the read length to {wildcards.readlen} bp."
    params:
        bptotrim=lambda wildcards: 120 - int(wildcards.readlen),
    shell:
        "(seqtk trimfq -e {params.bptotrim} {input} > "
        "raw_samples/seq_trim_{wildcards.species}sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq;"
        "seqtk seq -a raw_samples/seq_trim_{wildcards.species}sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq > "
        "{output}; "
        "rm raw_samples/seq_trim_{wildcards.species}sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l.fastq) 2> {log}"


rule deamsim:
    input:
        "raw_samples/seq_trim_{species}sp_{readlen}bp_{damage}d_{overhang}l.fasta",
    log:
        "raw_samples/deamsim_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        temp("raw_samples/deamsim_{species}sp_{readlen}bp_{damage}d_{overhang}l.fasta"),
    conda:
        "../envs/gargammel.yaml"
    message:
        "Changing the deamination rate to {wildcards.damage} and the overhang length to {wildcards.overhang}."
    params:
        nick=lambda wildcards: 0.03 if float(wildcards.damage) != 0 else 0,
        over=lambda wildcards: 0.4 if float(wildcards.damage) != 0 else 0,
        deam2=lambda wildcards: 0.01 if float(wildcards.damage) != 0 else 0
    shell:
        "(deamSim -damage {params.nick},{params.over},{params.deam2},{wildcards.damage} {input} > {output}) 2> {log}"



rule art_illumina:
    input:
        "raw_samples/deamsim_{species}sp_{readlen}bp_{damage}d_{overhang}l.fasta",
    log:
        "raw_samples/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        "raw_samples/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.fastq.gz",
    conda:
        "../envs/gargammel.yaml"
    message:
        "Simulating a library for sample effect_{wildcards.species}sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l."
    params:
        art_len=lambda wildcards: int(wildcards.readlen) if int(wildcards.readlen) != 120 else 100,
        out_basename = "raw_samples/illumina_{species}sp_{readlen}bp_{damage}d_{overhang}l"
    shell:
        "(art_illumina -ss HS25 -amp -na -i {input} -l {params.art_len} -c 1 -qs 0 -qs2 0 "
        "-o {params.out_basename}; bgzip --force --stdout {params.out_basename}.fq > {output}; "
        "rm {params.out_basename}.fq) 2> {log}"



