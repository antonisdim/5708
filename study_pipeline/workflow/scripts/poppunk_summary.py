#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import os


def poppunk_summary(input_file, out_primary_clusters, out_cluster_counts):
    """Function to get the primary cluster assignment for each assembly and
    count the number of individuals assigned to each cluster."""

    # read the input file
    poppunk_dump = pd.read_csv(input_file, sep=",")

    # cleanup the sample name
    poppunk_dump["Taxon"] = (
        poppunk_dump["Taxon"]
        .apply(lambda x: os.path.basename(x))
        .str.replace("_scaffolds.fasta", "")
    )

    # give unique ids to each cluster
    poppunk_dump["Cluster_id"] = poppunk_dump.groupby("Cluster").ngroup().add(1)

    # get the sample and its primary cluster
    primary_clusters = poppunk_dump[["Taxon", "Cluster_id"]].rename(
        columns={"Cluster_id": "Cluster"}
    )

    # store primary clusters to file
    primary_clusters.to_csv(out_primary_clusters, sep="\t", header=True, index=False)

    # groupby cluster, counts the number of assigned individuals, and write to file
    primary_clusters.groupby("Cluster").count().reset_index().rename(
        columns={"Taxon": "Count"}
    ).sort_values(by=["Count"], ascending=False).to_csv(
        out_cluster_counts, sep="\t", index=False, header=True
    )


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    poppunk_summary(
        input_file=snakemake.input[0],
        out_primary_clusters=snakemake.output.cluster_ids,
        out_cluster_counts=snakemake.output.cluster_counts,
    )
