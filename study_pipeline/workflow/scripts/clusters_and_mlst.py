#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from utilities import get_mlst_header


def clusters_and_mlst(clusters, mlst, out_file):
    """Function to count how many STs in each cluster"""

    # get correct header for mlst
    mlst_header = get_mlst_header(mlst)

    # read inputs
    cluster_table = pd.read_csv(clusters, sep=",", names=["sample", "cluster"])
    mlst_table = pd.read_csv(
        mlst, sep="\t", names=mlst_header, usecols=["Sample", "ST"]
    )

    clust_vs_mlst = pd.merge(
        cluster_table,
        mlst_table.rename(columns={"Sample": "sample"}),
        how="inner",
        on="sample",
    )

    # group and write to file
    clust_vs_mlst.groupby(["cluster", "ST"]).count().to_csv(
        out_file, sep="\t", index=True, header=True
    )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    clusters_and_mlst(
        clusters=snakemake.input.clusters,
        mlst=snakemake.input.mlst,
        out_file=snakemake.output[0],
    )
