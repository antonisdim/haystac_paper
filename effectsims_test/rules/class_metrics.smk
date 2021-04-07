#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule extract_read_names:
    input:
        get_tax_frags,
    log:
        "class_analysis_{mis}_{thresh}/ROC_{species}taxa_reads.log",
    output:
        "class_analysis_{mis}_{thresh}/ROC_{species}taxa_reads.tsv",
    message:
        "Getting a list of all the reads and the species they come from that were used for the gargammel sims."
    script:
        "../scripts/extract_read_names.py"


rule roc_per_read:
    input:
        lkhd_matrix="effect_analysis_out_{mis}/probabilities/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l_likelihood_ts_tv_matrix.csv",
        taxonomy="class_analysis_{mis}/ROC_{species}taxa_reads.tsv",
    log:
        "class_analysis_{mis}_{thresh}/Class_metrics_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        "class_analysis_{mis}_{thresh}/Class_metrics_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.tsv",
    message:
        "Calculating precision, recall and f1 score for matrix {input}."
    resources:
        concurrent_roc=1,
    script:
        "../scripts/class_metrics.py"


rule grid_roc_per_read:
    input:
        lkhd_matrix="effect_grid_analysis_out_{mis}_{thresh}/probabilities/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l/effect_{species}sp_{readlen}bp_{damage}d_{overhang}l_likelihood_ts_tv_matrix.csv",
        taxonomy="class_analysis_{mis}_{thresh}/ROC_{species}taxa_reads.tsv",
    log:
        "class_analysis_{mis}_{thresh}/Class_metrics_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.log",
    output:
        "class_analysis_{mis}_{thresh}/Class_metrics_effect_{species}sp_{readlen}bp_{damage}d_{overhang}l.tsv",
    message:
        "Calculating precision, recall and f1 score for matrix {input}."
    resources:
        concurrent_roc=1,
    script:
        "../scripts/class_metrics.py"
