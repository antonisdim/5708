#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import csv

SAMPLE_TABLE = "samples.tsv"
REF_GENOME_TABLE = "reference_genomes.tsv"


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
    if hasattr(wildcards, "cluster"):
        fastbaps = checkpoints.run_fastbaps.get(pathogen=wildcards.pathogen)
        clusters = pd.read_csv(
            fastbaps.output[0], sep=",", names=["Sample_Acc", "Cluster"]
        )
        cluster_samples = clusters[(clusters["Cluster"] == int(wildcards.cluster))]
        patho_samples = patho_samples[
            patho_samples["Sample_Acc"].isin(cluster_samples["Sample_Acc"])
        ]

    inputs_all = []

    for key, sam in patho_samples.iterrows():
        inputs_all.append(sam["Sample_Acc"])

    return inputs_all


def get_mlst_header(input_file):
    """Function to get the correct mlst header"""

    # get the number of columns
    with open(input_file) as fin:
        reader = csv.reader(fin, delimiter="\t")
        first_row = next(reader)
        num_cols = len(first_row)

    # start the header
    st_info = ["Sample", "Scheme", "ST"]
    loci = []

    # fill in the number of mlst loci
    if len(st_info) < num_cols:
        loci = [f"Locus_{num + 1}" for num in range(num_cols - len(st_info))]

    # full header
    header = st_info + loci

    return header


def get_ref_genome(wildcards):
    """Function to get the right ref genome for each fastbaps cluster"""

    ref_genomes = pd.read_csv(
        REF_GENOME_TABLE,
        sep="\t",
    )

    return ref_genomes.loc[
        ref_genomes["Cluster"] == int(wildcards.cluster), "Accession"
    ][0]