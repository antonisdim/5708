#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

SAMPLE_TABLE = "samples.tsv"


def read_sample_list():
    """Read the user sample table"""

    samples = pd.read_csv(
        SAMPLE_TABLE,
        sep="\t",
        names=["Sample_Acc", "Species", "Data_source"],
    )

    return samples


def get_right_pathogen(wildcards, checkpoints):
    """Get the E. coli samples"""

    samples = read_sample_list()

    species = ""
    if wildcards.pathogen == "ecoli":
        species = "Escherichia coli"

    patho_samples = samples[(samples["Species"] == species)]

    # check if we need context or not
    if hasattr(wildcards, "dataset"):
        if wildcards.dataset == "no":
            patho_samples = samples[
                (samples["Species"] == species) & (samples["Data_source"] == "SB27")
            ]
        else:
            patho_samples = samples[(samples["Species"] == species)]

    # check what cluster if necessary
    if hasattr(wildcards, "num"):
        poppunk = checkpoints.process_poppunk.get(
            pathogen=wildcards.pathogen, dataset=wildcards.dataset
        )
        primary_clusters = pd.read_csv(poppunk.output.cluster_ids, sep="\t")
        cluster_samples = primary_clusters[
            (primary_clusters["Cluster"] == int(wildcards.num))
        ]
        patho_samples = patho_samples[
            patho_samples["Sample_Acc"].isin(cluster_samples["Taxon"])
        ]

    inputs_all = []

    for key, sam in patho_samples.iterrows():
        inputs_all.append(sam["Sample_Acc"])

    return inputs_all
