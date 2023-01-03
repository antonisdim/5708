#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import csv
import os

SAMPLE_TABLE = "aux_files/samples.tsv"
CONT_SAMPLE_TABLE = "aux_files/excluded_samples.tsv"
LINEAGE_CONTX_TABLE = "aux_files/lineage_context_samples.tsv"
REF_GENOME_TABLE = "aux_files/reference_genomes.tsv"
OUT_GENOME_TABLE = "aux_files/outgroup_genomes.tsv"


def read_sample_list():
    """Read the user sample table"""

    samples = pd.read_csv(
        SAMPLE_TABLE,
        sep="\t",
        names=["Sample_Acc", "Species", "Data_source"],
    )

    return samples


def read_lineage_contx_list():
    """Read the user lineage context sample table"""

    samples = pd.read_csv(
        LINEAGE_CONTX_TABLE,
        sep="\t",
        names=["Sample_Acc", "Species", "Data_source"],
        dtype=object,
    )

    return samples


def get_contaminated_samples():
    """Function to get the contaminated samples in our dataset"""

    contaminated = pd.DataFrame(columns=["Sample_Acc", "Species"])

    if CONT_SAMPLE_TABLE:
        contaminated = pd.read_csv(
            CONT_SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species"]
        )

    return contaminated


def get_right_pathogen(wildcards, checkpoints, cont_check=True):
    """Get the E. coli samples"""

    samples = read_sample_list()

    species = ""
    genus = ""

    if wildcards.pathogen == "ecoli":
        species = "Escherichia coli"
    if wildcards.pathogen == "campylobacter":
        genus = "Campylobacter"

    if species:
        patho_samples = samples[(samples["Species"] == species)]
    else:
        patho_samples = samples[(samples["Species"].str.contains(genus, na=False))]

    # check if we need context or not
    if hasattr(wildcards, "dataset"):
        if wildcards.dataset == "no":
            patho_samples = patho_samples[(patho_samples["Data_source"] == "SB27")]

    # check what cluster if necessary
    if hasattr(wildcards, "cluster") and (wildcards.cluster != "1000"):
        fastbaps = checkpoints.run_fastbaps.get(pathogen=wildcards.pathogen)
        clusters = pd.read_csv(
            fastbaps.output[0], sep=",", names=["Sample_Acc", "Cluster"]
        )
        cluster_samples = clusters[(clusters["Cluster"] == int(wildcards.cluster))]
        patho_samples = patho_samples[
            patho_samples["Sample_Acc"].isin(cluster_samples["Sample_Acc"])
        ]

    # check if I need to add global context to a specific lineage
    if hasattr(wildcards, "population") and (wildcards.population == "population"):
        lineage_context = read_lineage_contx_list()
        lineage_context = lineage_context[
            (lineage_context["Data_source"] == wildcards.cluster)
        ]
        patho_samples = pd.concat([patho_samples, lineage_context], ignore_index=True)

    # exclude any contaminated samples
    if cont_check:
        contaminated_samples = get_contaminated_samples()
        patho_samples = patho_samples[
            ~patho_samples["Sample_Acc"].isin(contaminated_samples["Sample_Acc"])
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
        (ref_genomes["Cluster"] == int(wildcards.cluster))
        & (ref_genomes["Species"] == wildcards.pathogen),
        "Accession",
    ].to_list()[0]


def genome_chromosome(taxon_fasta_idx):
    """Get the accession for the main chromosome"""

    faidx = pd.read_csv(
        taxon_fasta_idx,
        sep="\t",
        names=["Name", "Length", "Offset", "Linebases", "Linewidth"],
    )

    chromosome = faidx.sort_values(by="Length", ascending=False).iloc[0, 0]

    return chromosome


def get_out_genome(wildcards):
    """Function to get the right ref genome for each fastbaps cluster"""

    ref_genomes = pd.read_csv(
        OUT_GENOME_TABLE,
        sep="\t",
    )

    return ref_genomes.loc[
        (ref_genomes["Cluster"] == int(wildcards.cluster))
        & (ref_genomes["Species"] == wildcards.pathogen),
        "Accession",
    ].to_list()[0]


def get_read_file(wildcards, mate):
    """Get the relative path to the R1 read file"""

    sb27_samples = read_sample_list()

    if wildcards.sample in sb27_samples["Sample_Acc"].values:
        sam_path = f"adRm/{wildcards.sample}_R{mate}_adRm.fastq.gz"
    else:
        sam_path = f"sra_context_ecoli/{wildcards.sample}_R{mate}_adRm.fastq.gz"

    return sam_path


def get_r1(wildcards):
    """Get the relative path to the R1 read file"""

    return get_read_file(wildcards, "1")


def get_r2(wildcards):
    """Get the relative path to the R2 read file"""

    return get_read_file(wildcards, "2")


def get_clusters_to_run(wildcards):
    """Get the clusters that summary stats will be run on"""

    ref_genomes = pd.read_csv(
        REF_GENOME_TABLE,
        sep="\t",
    )
    ref_genomes = ref_genomes.loc[
        ref_genomes["Species"] == wildcards.pathogen,
    ]

    clusters = ref_genomes["Cluster"].unique().tolist()
    clusters.pop(0)
    clusters.pop(-1)  # todo remove that

    return clusters


def get_ref_idx(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{wildcards.pathogen}/{ref}.fasta.gz.fai"


def population_host_metadata(pop_metadata):
    # read the sample metadata file
    pop_meta = pd.read_csv(pop_metadata, sep="\t")

    # define lists with the broader host category
    human = ["Human"]
    ruminant = ["Beef", "Bovine", "Sheep", "Goat", "Ovine/Goat"]
    avian = ["Chicken", "Avian", "Crow", "Poultry", "Turkey"]
    food = ["Food", "Dairy"]
    swine = ["Pork", "Swine"]
    other_mammal = ["Primate", "Rodent", "Deer", "Canine"]
    other = [
        "Soil",
        "Water",
        "Water/River",
        "Soil/Dust",
        "ND/Other",
        "Laboratory",
        "Plant",
        "Animal-related",
    ]

    # assign trait value based on host organism/type
    pop_meta.loc[pop_meta["Host"].isin(human), "Trait"] = "Human"
    pop_meta.loc[pop_meta["Host"].isin(ruminant), "Trait"] = "Ruminant"
    pop_meta.loc[pop_meta["Host"].isin(avian), "Trait"] = "Avian"
    pop_meta.loc[pop_meta["Host"].isin(food), "Trait"] = "Food"
    pop_meta.loc[pop_meta["Host"].isin(swine), "Trait"] = "Swine"
    pop_meta.loc[pop_meta["Host"].isin(other_mammal), "Trait"] = "Other_mammal"
    pop_meta.loc[pop_meta["Host"].isin(other), "Trait"] = "Other"

    return pop_meta
