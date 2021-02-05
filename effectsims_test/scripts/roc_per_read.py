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
    """Function to calculate the roc"""

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

    aln_reads_vector = (
        ts_tv_matrix[["Taxon", "Read_ID"]]
        .groupby("Taxon")
        .count()
        .squeeze()
        .rename("Taxon")
    )
    aln_reads_vector["Dark_Matter"] = 0
    aln_reads_vector["Grey_Matter"] = 0

    # Sum the Dirichlet Assignments per taxon and calculate the Dark Matter reads from the Dirichlet Assignment column
    ts_tv_group = ts_tv_matrix.groupby("Read_ID").sum().squeeze()
    grey_matter = ts_tv_group.where(ts_tv_group == 0).replace(0, 1).fillna(0)

    if len(ts_tv_matrix.Taxon.unique()) > 1:
        a = ts_tv_matrix.groupby("Taxon").sum().squeeze().astype(float)
    else:
        a = ts_tv_matrix.groupby("Taxon").sum().iloc[:, 0].astype(float)
    a.loc["Grey_Matter"] = grey_matter.sum()

    # Add the non aligned filtered reads count in the Dark Matter category
    total_fastq_reads = len(taxa["reads"].unique())
    reads_in_bams = len(ts_tv_matrix["Read_ID"].unique())

    remaining_dark_matter = total_fastq_reads - reads_in_bams

    a.loc["Dark_Matter"] = remaining_dark_matter

    # split them in assigned and unassigned reads
    assigned_reads = ts_tv_matrix[ts_tv_matrix["Dirichlet_Assignment"] != 0]
    assigned_comp = pd.merge(
        assigned_reads, taxa, left_on=["Read_ID"], right_on=["reads"]
    )

    # unassigned_reads = ts_tv_matrix[ts_tv_matrix["Dirichlet_Assignment"] == 0]
    # unassigned_comp = pd.merge(
    #     unassigned_reads, taxa, left_on=["Read_ID"], right_on=["reads"]
    # )
    #
    # unassigned_comp["TN"] = np.where(
    #     (unassigned_comp["Taxon"] != unassigned_comp["species"]), 1, 0
    # )

    # start the counter for false negative, false positive, true positive and true negative
    fn_count = a.loc["Dark_Matter"] + a.loc["Grey_Matter"]
    fp_count = 0
    tp_count = 0
    # tn_count = unassigned_comp["TN"].sum()
    # true negative reads do not exist as all of the reads have a species in the database
    tn_count = 0

    # merge the tables

    db_species = assigned_reads["Taxon"].unique().tolist()

    for index, row in assigned_comp.iterrows():
        if row["Taxon"] == row["species"]:
            tp_count += 1
        else:
            # if simulated species not in the database
            if row["species"] not in db_species:
                fp_count += 1
            # if simulated species in the database
            else:
                fp_count += 1
                fn_count += 1

    # tp_rate = tp_count / (tp_count + fn_count)
    # fp_rate = fp_count / (fp_count + tn_count)
    # fn_rate = fn_count / (fn_count + tp_count)
    #
    # # print the rates to the output
    # with open(output, "w") as fout:
    #     print(tp_rate, fp_rate, fn_rate, sep="\t", file=fout)

    with open(output, "w") as fout:
        print(tp_count, fp_count, fn_count, sep="\t", file=fout)

    print(
        "Done with the rates.", file=sys.stderr,
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
