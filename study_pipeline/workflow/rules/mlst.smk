#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_right_pathogen


rule mlst:
    input:
        "assemblies/{accession}_scaffolds.fasta",
    log:
        "mlst/{accession}_SB27_{pathogen}_{dataset}_contx.log",
    output:
        temp("mlst/{accession}_SB27_{pathogen}_{dataset}_contx.tsv"),
    message:
        "Running MLST typing for sample {wildcards.accession} for {wildcards.pathogen}, "
        "{wildcards.dataset} context."
    conda:
        "../envs/mlst.yaml"
    wildcard_constraints:
        dataset="(no|plus)",
    shell:
        "mlst {input} --scheme {wildcards.pathogen} --nopath --label {wildcards.accession} 1> {output} 2> {log}"


def get_mlst_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(
            f"mlst/{sample}_SB27_{wildcards.pathogen}_{wildcards.dataset}_contx.tsv"
        )

    return input_paths


rule aggregate_mlst:
    input:
        get_mlst_paths,
    output:
        "mlst/SB27_{pathogen}_{dataset}_contx_cluster_{num}.tsv",
    message:
        "Aggregating all the MLST results for {wildcards.pathogen} from SB27, with {wildcards.dataset} "
        "context, for PopPUNK cluster {wildcards.num}."
    shell:
        "cat {input} > {output}"


rule count_mlst:
    input:
        "mlst/SB27_{pathogen}_{dataset}_contx_cluster_{num}.tsv",
    output:
        "mlst/SB27_{pathogen}_{dataset}_contx_cluster_{num}_mlst_count.tsv",
    message:
        "Counting the number of individuals for cluster {wildcards.num} of {wildcards.pathogen}, from "
        "SB27, {wildcards.dataset} context."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count_mlst.py"
