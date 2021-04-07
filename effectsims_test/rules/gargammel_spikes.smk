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


rule gargammel_spike_inputs:
    input:
        "spike_accs/db_taxa_accessions.tsv",
    output:
        "garg_sims_spikes/{taxon}/endo/{accession}.fasta",
    message:
        "Preparing the gargammel inputs for species {wildcards.taxon} with accession {wildcards.accession}."
    shell:
        "mkdir -p garg_sims_spikes/{wildcards.taxon}/endo garg_sims_spikes/{wildcards.taxon}/bact garg_sims_spikes/{wildcards.taxon}/cont;"
        "gzip -c -d {CACHE}/ncbi/{wildcards.taxon}/{wildcards.accession}.fasta.gz > {output}"


def get_tax_acc_spike_pairs(_):
    """Function to get the pairs of tax and acc for the gargammel spiked sims"""

    spike_scc = checkpoints.download_spike_accs.get()
    tax_acc = pd.read_csv(
        spike_scc.output.sim_list, sep="\t", names=["Taxon", "Accession"]
    )

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row["Taxon"]
        accession = row["Accession"]
        paths.append(f"garg_sims_spikes/{taxon}/endo/{accession}.fasta")

    return paths


rule all_spike_dirs_prepared:
    input:
        get_tax_acc_spike_pairs,
    output:
        "garg_sims_spikes/all_spike_dirs_prepared.done",
    message:
        "All gargammel inputs have been prepared."
    shell:
        "touch {output}"


rule simulate_gargammel_spike:
    input:
        "garg_sims_spikes/all_spike_dirs_prepared.done",
        misince="garg_sim_param_files/ERR2181127_dnacomp.txt",
        frag_len="garg_sim_param_files/ERR2181127_lgdistribution.txt",
        mapdamagee="garg_sim_param_files/ERR2181127_misincorporation.txt",
    log:
        "garg_sims_spikes/{taxon}.log",
    output:
        "garg_sims_spikes/{taxon}_s.fq.gz",
    message:
        "Simulating 100 reads from taxon {wildcards.taxon}."
    conda:
        "../envs/gargammel.yaml"
    shell:
        "( gargammel --comp 0,0,1 -n 100 --misince {input.misince} -f {input.frag_len} --mapdamagee {input.mapdamagee} "
        "single -rl 125 -se -ss HS25 garg_sims_spikes/{wildcards.taxon}/ -o garg_sims_spikes/{wildcards.taxon} ) &> {log}"


rule create_reaL_human_spiked_sim:
    input:
        spike="garg_sims_spikes/{taxon}_s.fq.gz",
        human="real_human/ERR2181127.fastq.gz"
    output:
        temp("raw_samples/real_spiked_lib_{taxon}.fastq.gz"),
    message:
        "Spiking library ERR2181127 with 100 reads from {wildcards.taxon}.",
    shell:
        "cat {input.spike} {input.human} 1> {output}"


rule adapterremoval_single_end:
    input:
        fastq="raw_samples/real_spiked_lib_{taxon}.fastq.gz",
    log:
        "adRm_spikes/real_spiked_lib_{taxon}_adRm.log",
    output:
        "adRm_spikes/real_spiked_lib_{taxon}_adRm.fastq.gz",
    message:
        "Trimming sequencing adapters from file {input.fastq}."
    conda:
        "../envs/adapterremoval.yaml"
    threads: 5
    params:
        basename="adRm_spikes/real_spiked_lib_{taxon}_adRm",
    shell:
        "(AdapterRemoval"
        "   --file1 {input}"
        "   --basename {params.basename} "
        "   --gzip "
        "   --minlength 15 "
        "   --threads {threads}"
        "   --trimns && "
        " mv {params.basename}.truncated.gz {output}"
        ") 2> {log}"


def get_real_human_lib_paths(_):
    """Function to get the libraries for the human spiked sims"""

    spike_scc = checkpoints.download_spike_accs.get()
    tax_acc = pd.read_csv(
        spike_scc.output.sim_list, sep="\t", names=["Taxon", "Accession"]
    )

    paths = []

    for index, row in tax_acc.iterrows():
        taxon = row["Taxon"]
        paths.append(f"adRm_spikes/real_spiked_lib_{taxon}_adRm.fastq.gz")

    return paths


rule all_real_human_spiked:
    input:
        get_real_human_lib_paths,
    output:
        "all_real_human_spiked.done",
    shell:
        "touch {output}"
