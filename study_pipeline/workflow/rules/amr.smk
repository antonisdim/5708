#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

os.makedirs("ariba_tmp", exist_ok=True)


rule ariba_getref:
    log:
        "public_dbs/get_{db}.log",
    output:
        "public_dbs/{db}.fa",
        "public_dbs/{db}.tsv",
    message:
        "Getting the fasta file for the {wildcards.db} database."
    conda:
        "../envs/amr.yaml"
    params:
        out_basename="public_dbs/{db}",
    shell:
        "ariba getref card {params.out_basename} &> {log}"


rule ariba_prepareref:
    input:
        fasta="public_dbs/{db}.fa",
        meta="public_dbs/{db}.tsv",
    output:
        "ariba_{db}/00.info.txt",
        "ariba_{db}/01.filter.check_metadata.tsv",
        "ariba_{db}/02.cdhit.clusters.tsv",
    message:
        "Preparing the {wildcards.db} database for ariba."
    conda:
        "../envs/amr.yaml"
    threads: workflow.cores
    params:
        cdhit_mem=0,
        outdir="ariba_{db}",
    shell:
        "ariba prepareref --fasta {input.fasta} --metadata {input.meta} --threads {threads} "
        "--cdhit_max_memory {params.cdhit_mem} --force {params.outdir}"


rule ariba_run:
    input:
        r1="adRm/{sample}_R1_adRm.fastq.gz",
        r2="adRm/{sample}_R2_adRm.fastq.gz",
        db="ariba_{db}/02.cdhit.clusters.tsv",
    output:
        "ariba_{db}_{pathogen}/{sample}/assembled_genes.fa.gz",
        "ariba_{db}_{pathogen}/{sample}/assembled_seqs.fa.gz",
        "ariba_{db}_{pathogen}/{sample}/assemblies.fa.gz",
        "ariba_{db}_{pathogen}/{sample}/debug.report.tsv",
        "ariba_{db}_{pathogen}/{sample}/log.clusters.gz",
        "ariba_{db}_{pathogen}/{sample}/report.tsv",
    message:
        "Running ariba on sample {wildcards.sample} against the {wildcards.db} db."
    conda:
        "../envs/amr.yaml"
    threads: 1
    params:
        db="ariba_{db}",
        outdir="ariba_{db}_{pathogen}/{sample}",
        assembler="spades",
        tmpdir="ariba_tmp",
    shell:
        "ariba run {params.db} {input.r1} {input.r2} {params.outdir} --assembler {params.assembler} "
        "--tmp_dir {params.tmpdir} --threads {threads} --force"


def get_cluster_reports(wildcards):
    """Get the samples belonging to a cluster to summarise the ariba reports"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(
            f"ariba_{wildcards.db}_{wildcards.pathogen}/{sample}/report.tsv"
        )

    return input_paths


rule ariba_summary:
    input:
        get_cluster_reports,
    output:
        "ariba_{db}_{pathogen}/{pathogen}_cluster_{cluster}_report.csv",
    message:
        "Summarising reports for {wildcards.pathogen} cluster {wildcards.cluster} for {wildcards.db} db."
    conda:
        "../envs/amr.yaml"
    params:
        out_prefix="ariba_{db}_{pathogen}/{pathogen}_cluster_{cluster}_report",
    shell:
        "ariba summary {params.out_prefix} {input} --preset all"
