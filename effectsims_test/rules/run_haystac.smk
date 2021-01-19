#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule run_haystac_sample:
    input:
        "raw_samples/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.fastq.gz",
    log:
        "haystac_sample/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.log"
    output:
        "haystac_sample/effect_5652sp_{readlen}bp_{damage}d_{overhang}l/fastq_inputs/meta/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.size"
    message:
        "Pre processing sample {input}."
    params:
        outdir_basename="effect_5652sp_{readlen}bp_{damage}d_{overhang}l"
    threads: 1
    shell:
        "haystac sample --fastq {input} --trim-adapters False "
        "--output haystac_sample/{params.outdir_basename} --cores {threads} &> {log}"


rule run_haystac_analyse:
    input:
        "haystac_sample/effect_5652sp_{readlen}bp_{damage}d_{overhang}l/fastq_inputs/meta/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.size",
    log:
        "effect_analysis_out/effect_5652sp_{readlen}bp_{damage}d_{overhang}l.log"
    output:
        "effect_analysis_out/probabilities/effect_5652sp_{readlen}bp_{damage}d_{overhang}l/effect_5652sp_{readlen}bp_{damage}d_{overhang}l_posterior_abundance.tsv"
    message:
        "Analysing sample effect_5652sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l."
    params:
        outdir_basename="effect_analysis_out"
    threads: 15
    resources:
        concurrent_samples=1
    shell:
        "haystac analyse --mode abundances --database paper_db "
        "--sample haystac_sample/effect_5652sp_{wildcards.readlen}bp_{wildcards.damage}d_{wildcards.overhang}l "
        "--cores {threads} --output {params.outdir_basename} &> {log}"