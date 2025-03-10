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
        'scheme=$"{wildcards.pathogen}";'
        'if [[ "{wildcards.pathogen}" == "cjejuni" ]]; then'
        '   scheme=$"campylobacter";'
        "   mlst {input} --scheme $scheme --nopath --label {wildcards.accession} 1> {output} 2> {log};"
        "else"
        "   mlst {input} --scheme {wildcards.pathogen} --nopath --label {wildcards.accession} 1> {output} 2> {log};"
        "fi"


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
        "mlst/SB27_{pathogen}_{dataset}_contx.tsv",
    message:
        "Aggregating all the MLST results for {wildcards.pathogen} from SB27, with {wildcards.dataset} "
        "context."
    shell:
        "cat {input} > {output}"


rule count_mlst:
    input:
        "mlst/SB27_{pathogen}_{dataset}_contx.tsv",
    output:
        "mlst/mlst_summary_SB27_{pathogen}_{dataset}_contx.tsv",
    message:
        "Counting the number of individuals of {wildcards.pathogen} for each ST, from SB27, "
        "{wildcards.dataset} context."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count_mlst.py"
