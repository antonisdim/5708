#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import os


def poppunk_summary(input_file, out_cluster_counts):
    """Function to get the primary cluster assignment for each assembly and
    count the number of individuals assigned to each cluster."""

    # read the input file
    fastbaps_clusters = pd.read_csv(input_file, sep=",", names=["sample", "cluster"])

    # groupby cluster, counts the number of assigned individuals
    cluster_counts = (
        fastbaps_clusters.groupby("cluster")
        .count()
        .reset_index()
        .sort_values(by="sample", ascending=False)
    )

    # percentage of assigned samples to each cluster
    cluster_counts["percent"] = (
        cluster_counts["sample"] / cluster_counts["sample"].sum()
    ) * 100

    # make it clean
    clean_counts = cluster_counts.round(2).rename(
        columns={"cluster": "Cluster", "sample": "Sample Count", "percent": "Sample %"}
    )
    # write to file
    clean_counts.to_csv(out_cluster_counts, sep="\t", header=True, index=False)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    poppunk_summary(
        input_file=snakemake.input[0],
        out_cluster_counts=snakemake.output[0],
    )
