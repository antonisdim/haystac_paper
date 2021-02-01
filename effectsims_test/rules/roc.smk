#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule extract_acc_fasta:
    input:
        get_tax_acc_pairs,
    log:
        "roc_analysis_{mis}/ROC_taxa_acc.log",
    output:
        "roc_analysis_{mis}/ROC_taxa_acc.tsv",
    message:
        "Getting a list of all the accessions that exist in every fasta file used for the gargammel sims."
    script:
        "../scripts/extract_acc_fasta.py"


rule roc_per_read:
    input:
        lkhd_matrix="effect_analysis_out_{mis}/probabilities/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l_likelihood_ts_tv_matrix.csv",
        taxonomy="roc_analysis/ROC_taxa_acc.tsv",
    log:
        "roc_analysis_{mis}/ROC_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        "roc_analysis_{mis}/ROC_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.tsv",
    message:
        "Calculating the ROC for matrix {input}."
    script:
        "../scripts/roc_per_read.py"
