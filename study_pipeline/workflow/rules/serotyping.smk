#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_right_pathogen, get_ref_genome


rule serotyping:
    input:
        "seqs_{pathogen}/{sample}_ref_{accession}.fasta",
    log:
        "serotyping_{pathogen}/{sample}_ref_{accession}.log",
    output:
        table=temp("serotyping_{pathogen}/{sample}_ref_{accession}.tsv"),
        dir=temp(directory("serotyping_{pathogen}/{sample}_ref_{accession}")),
    message:
        "Getting the surface antigen types for {wildcards.pathogen} sample {wildcards.sample}."
    conda:
        "../envs/ectyper.yaml"
    shell:
        "(ectyper --input {input} --output {output.dir} && mv {output.dir}/output.tsv {output.table}) &> {log}"


def get_serotyping_cluster_reports(wildcards):
    """Get the samples belonging to a cluster to summarise the O and H surface antigen reports"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(f"serotyping_{wildcards.pathogen}/{sample}_ref_{ref}.tsv")

    return input_paths


rule serotyping_summary:
    input:
        get_serotyping_cluster_reports,
    output:
        "serotyping_{pathogen}/{pathogen}_{population}_{cluster}_sero_report.tsv",
    message:
        "Summarising the surface antigen reports for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/ectyper.yaml"
    params:
        ref=lambda wildcards: get_ref_genome(wildcards),
    shell:
        "awk 'FNR>1 || NR==1' {input} | sed 's/_ref_{params.ref}//g' > {output}"
