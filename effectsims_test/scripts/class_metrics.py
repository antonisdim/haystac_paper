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


def roc_per_read(taxonomy, lkhd_matrix, output):
    """Function to calculate the classification metrics"""

    print(
        "Calculating the TP, FP and FN rates ...", file=sys.stderr,
    )

    # load the list of reads and the species they came from

    taxa = pd.read_csv(taxonomy, names=["species", "reads"], sep="\t")
    taxa["species"] = taxa["species"].str.replace(" ", "_")
    taxa["reads"] = taxa["reads"] + "-1"

    # load ts and tv matrix and calculate dark and grey matter

    ts_tv_matrix = pd.read_csv(
        lkhd_matrix, sep=",", usecols=["Taxon", "Read_ID", "Dirichlet_Assignment"]
    )

    # Small df for the grey matter reads
    ts_tv_group = ts_tv_matrix.groupby("Read_ID").sum().squeeze()
    grey_matter = ts_tv_group.where(ts_tv_group == 0).replace(0, 1).fillna(0)
    grey_df = grey_matter[grey_matter == 1].to_frame().reset_index()
    grey_df["Taxon"] = "Grey_Matter"

    # Small df for the dark matter reads
    dark_df = taxa[~taxa["reads"].isin(ts_tv_matrix["Read_ID"])].rename(
        columns={"reads": "Read_ID", "species": "Taxon"}
    )
    dark_df["Dirichlet_Assignment"] = 1
    dark_df["Taxon"] = "Dark_Matter"

    # split them in assigned and unassigned reads
    tax_assigned_reads = ts_tv_matrix[ts_tv_matrix["Dirichlet_Assignment"] != 0]
    assigned_reads = pd.concat([tax_assigned_reads, grey_df, dark_df], axis=0)

    # merge all the reads no matter where they were assigned
    assigned_comp = pd.merge(
        assigned_reads, taxa, left_on=["Read_ID"], right_on=["reads"]
    )

    # calculate classification stats

    values_true = assigned_comp["species"]
    values_pred = assigned_comp["Taxon"]

    precision = metrics.precision_score(
        values_true, values_pred, average="weighted", zero_division=0
    )
    recall = metrics.recall_score(
        values_true, values_pred, average="weighted", zero_division=0
    )
    f1 = metrics.f1_score(values_true, values_pred, average="weighted", zero_division=0)

    with open(output, "w") as fout:
        print(round(precision, 3), round(recall, 3), round(f1, 3), sep="\t", file=fout)

    print(
        "Done with the classification metrics.", file=sys.stderr,
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
