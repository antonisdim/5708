#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import csv

SAMPLE_TABLE = "aux_files/samples.tsv"
CONT_SAMPLE_TABLE = "aux_files/excluded_samples.tsv"
LINEAGE_CONTX_TABLE = "aux_files/lineage_context_samples.tsv"
REF_GENOME_TABLE = "aux_files/reference_genomes.tsv"
OUT_GENOME_TABLE = "aux_files/outgroup_genomes.tsv"
TIME_INTERVAL = 4
CLOCK_RATE_TABLE = "aux_files/clock_rates.tsv"


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


def get_contaminated_samples(wildcards):
    """Function to get the contaminated samples in our dataset"""

    contaminated = pd.DataFrame(columns=["Sample_Acc", "Species", "Pop_level"])

    if CONT_SAMPLE_TABLE:
        contaminated = pd.read_csv(
            CONT_SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species", "Pop_level"]
        )

        if hasattr(wildcards, "population") and (wildcards.population == "cluster"):
            contaminated = contaminated[
                contaminated["Pop_level"] == "Baps_Clusters"
            ].copy()

    return contaminated


def get_right_pathogen(wildcards, checkpoints, cont_check=True):
    """Get the E. coli samples"""

    samples = read_sample_list()

    species = ""
    genus = ""

    if wildcards.pathogen == "ecoli":
        species = "Escherichia coli"
    if wildcards.pathogen == "ccoli":
        species = "Campylobacter coli"
    if wildcards.pathogen == "cjejuni":
        species = "Campylobacter jejuni"
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
        contaminated_samples = get_contaminated_samples(wildcards)
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


def genome_chromosome(taxon_fasta_idx, size=False):
    """Get the accession for the main chromosome"""

    faidx = pd.read_csv(
        taxon_fasta_idx,
        sep="\t",
        names=["Name", "Length", "Offset", "Linebases", "Linewidth"],
    )

    if size:
        chromosome = faidx.sort_values(by="Length", ascending=False).iloc[0, 1]
    else:
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
    if len(clusters) > 1:
        clusters.pop(0)
    if wildcards.pathogen == "ecoli":  # todo remove that
        clusters.pop(-1)

    return clusters


def get_ref_idx(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{wildcards.pathogen}/{ref}.fasta.gz.fai"


def get_clock_rate(wildcards):
    """Get a clock rate estimation for a given taxon"""  # todo optimise this

    clock_rates = pd.read_csv(
        CLOCK_RATE_TABLE, sep="\t", names=["Species", "Cluster", "Rate"]
    )

    return clock_rates.loc[
        (clock_rates["Species"] == wildcards.pathogen)
        & (clock_rates["Cluster"] == int(wildcards.cluster)),
        "Rate",
    ].to_list()[0]


def population_host_metadata(pop_metadata):
    # read the sample metadata file
    pop_meta = pd.read_csv(pop_metadata, sep="\t")

    # define lists with the broader host category
    human = ["Human"]
    ruminant = ["Beef", "Bovine", "Sheep", "Goat", "Ovine/Goat"]
    avian = ["Chicken", "Avian", "Crow", "Poultry", "Turkey"]
    food = ["Food", "Dairy"]
    swine = ["Pork", "Swine", "Pig"]
    other_mammal = [
        "Primate",
        "Rodent",
        "Deer",
        "Canine",
        "Feline",
        "Equine",
        "Marine Mammal",
        "Other Mammal",
    ]
    other = [
        "Soil",
        "Water",
        "Water/River",
        "Soil/Dust",
        "ND/Other",
        "Laboratory",
        "Plant",
        "Animal-related",
        "Environment",
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


def get_aln_file(wildcards):
    """Get the right aln file for the corresponding summary statistic or tree"""

    whole_aln = (
        f"msa_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{wildcards.cluster}"
        f"_chr_aln.fasta"
    )
    whole_aln_nrec = (
        f"msa_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{wildcards.cluster}"
        f"_chr_aln_nrec.fasta"
    )
    snps_rec = (
        f"msa_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{wildcards.cluster}"
        f"_chr_aln_rec_snps.fasta"
    )
    snps_nrec = (
        f"msa_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{wildcards.cluster}"
        f"_chr_aln_nrec_snps.fasta"
    )

    # for the metrics
    if hasattr(wildcards, "metric"):
        if wildcards.metric in ["swfst", "tajima", "swfst-human", "tajima-human"]:
            input_path = whole_aln
        else:
            if wildcards.rec == "rec":
                input_path = snps_rec
            else:
                input_path = snps_nrec
    else:
        # for snp-sites at msa.smk
        if hasattr(wildcards, "rec"):
            if wildcards.rec == "rec":
                input_path = whole_aln
            else:
                input_path = whole_aln_nrec
        # for tree building
        else:
            if wildcards.cluster != "1000":
                input_path = snps_nrec
            else:
                input_path = snps_rec

    return input_path


def get_correct_metadata(wildcards):
    """Get the correct metadata file depending on context size"""

    if wildcards.population == "cluster":
        pop_meta = f"aux_files/{wildcards.pathogen}_all_meta.tsv"
    else:
        pop_meta = f"aux_files/{wildcards.pathogen}_lineage_all_meta.tsv"

    return pop_meta
