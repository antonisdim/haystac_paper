#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import numpy as np
import pandas as pd
import sys

from sklearn import metrics
from sys import argv


def roc_per_read(taxonomy, lkhd_matrix, output):
    """Function to calculate the roc"""

    print(
        "Calculating the ROC ...", file=sys.stderr,
    )

    taxa = pd.read_csv(taxonomy, names=["species", "acc"], sep="\t")
    taxa["species"] = taxa["species"].str.replace(" ", "_")

    # add human
    taxa = taxa.append(
        pd.DataFrame([["chr", "Homo_sapiens"]], columns=["acc", "species"])
    )

    lkhd = pd.read_csv(lkhd_matrix, sep=",")

    print("TSV file ready", file=sys.stderr)

    acc = "|".join(taxa.acc)
    lkhd["Accession"] = lkhd["Read_ID"].str.extract("(" + acc + ")", expand=False)
    lkhd = lkhd[~lkhd["Accession"].isnull()]
    complete_table = pd.merge(lkhd, taxa, left_on=["Accession"], right_on=["acc"])
    assigned_reads = complete_table
    assigned_reads["Truth"] = np.where(
        (assigned_reads["Taxon"] == assigned_reads["species"]), 1, 0
    )

    fpr, tpr, thresholds_roc = metrics.roc_curve(
        assigned_reads["Truth"], assigned_reads["Likelihood"]
    )
    roc = pd.DataFrame()
    roc["FPR"] = fpr
    roc["TPR"] = tpr
    roc["Thresholds"] = thresholds_roc
    roc.to_csv(output, sep="\t", index=None)

    print(
        "Done with the ROC.", file=sys.stderr,
    )


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    roc_per_read(
        taxonomy=snakemake.input.taxonomy,
        lkhd_matrix=snakemake.input.lkhd_matrix,
        output=snakemake.output[0],
    )
