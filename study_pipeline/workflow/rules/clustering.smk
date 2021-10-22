#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_right_pathogen


def get_sb27_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(f"assemblies/{sample}_scaffolds.fasta")

    return input_paths


rule create_poppunk_qfile:
    input:
        sample_list=get_sb27_paths,
    output:
        "poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_scaffold_paths.txt",
    message:
        "Creating the query file for PopPUNK for the SB27 {wildcards.pathogen} samples, "
        "{wildcards.dataset} context samples."
    wildcard_constraints:
        dataset="(no|plus)",
    shell:
        "for path in {input.sample_list}; do readlink -f $path; done 1> {output}"


rule run_poppunk:
    input:
        "poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_scaffold_paths.txt",
    log:
        "poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples/poppunk.log",
    output:
        "poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples/SB27_{dataset}_contx_{pathogen}_samples_clusters.csv",
    message:
        "Running PopPUNK for the SB27 samples of {wildcards.pathogen}, {wildcards.dataset} context samples."
    conda:
        "../envs/poppunk.yaml"
    threads: workflow.cores
    params:
        outdir="poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples",
        db="{pathogen}_poppunk",
    shell:
        "( poppunk --assign-query --ref-db {params.db} --q-files {input} --output {params.outdir} "
        "--ignore-length --threads {threads} ) 2> {log}"


checkpoint process_poppunk:
    input:
        "poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples/SB27_{dataset}_contx_{pathogen}_samples_clusters.csv",
    output:
        cluster_ids="poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples/"
        "SB27_{dataset}_contx_{pathogen}_samples_cluster_ids.tsv",
        cluster_counts="poppunk_{pathogen}/SB27_{dataset}_contx_{pathogen}_samples/"
        "SB27_{dataset}_contx_{pathogen}_samples_cluster_counts.tsv",
    message:
        "Quick summary of the PopPUNK cluster output for the SB27 samples of  {wildcards.pathogen}, "
        "{wildcards.dataset} context samples."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/poppunk_summary.py"
