#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from scripts.utilities import read_sample_list, get_right_pathogen


rule spades_assembly:
    input:
        fastq_r1="adRm/{accession}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{accession}_R2_adRm.fastq.gz",
    log:
        "assemblies/{accession}/spades.log",
    output:
        contigs="assemblies/{accession}_contigs.fasta",
        scaffolds="assemblies/{accession}_scaffolds.fasta",
    message:
        "Spades de novo assembly for {wildcards.accession}."
    conda:
        "../envs/spades.yaml"
    params:
        outdir="assemblies/{accession}",
    threads: 1
    shell:
        "( spades.py -1 {input.fastq_r1} -2 {input.fastq_r2} -o {params.outdir} --threads 1 --careful && "
        "mv {params.outdir}/contigs.fasta {output.contigs}; "
        "mv {params.outdir}/scaffolds.fasta {output.scaffolds} ) 2> {log}"


def get_assemblies(_):
    """Get the paths to the contigs"""

    samples = read_sample_list()[["Sample_Acc"]]
    inputs = []

    for key, sam in samples.iterrows():
        assembly = f"assemblies/{sam['Sample_Acc']}_contigs.fasta"
        inputs.append(assembly)

    return inputs


rule get_assemblies:
    input:
        get_assemblies,
    output:
        "assemblies.done",
    message:
        "Spades has finished performing the de novo assemblies."
    shell:
        "touch {output}"


rule run_assembly_stats:
    input:
        assembly="assemblies/{accession}_scaffolds.fasta",
    log:
        "assemblies/{accession}_scaffolds_stats.log",
    output:
        assembly_stats=temp("assemblies/{accession}_scaffolds_stats.tsv"),
    message:
        "Calculating assembly stats for sample {wildcards.accession}."
    params:
        extra="-t",
    threads: 1
    wrapper:
        "0.80.2/bio/assembly-stats"


def get_assemblies_stats(wildcards):
    """Get the paths to the contigs"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    inputs = []

    for sample in samples_list:
        assembly = f"assemblies/{sample}_scaffolds_stats.tsv"
        inputs.append(assembly)

    return inputs


rule get_assemblies_stats:
    input:
        get_assemblies_stats,
    output:
        "assemblies/assemblies_stats_{pathogen}.tsv",
    message:
        "Assembly stats have been calculated for {wildcards.pathogen}."
    shell:
        "awk 'FNR>1 || NR==1' {input} | sed 's/assemblies\///g' | sed 's/_scaffolds.fasta//g' > {output}"
