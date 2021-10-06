#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

SAMPLE_TABLE = "samples.tsv"


rule spades_assembly:
    input:
        fastq_r1="adRm/{accession}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{accession}_R2_adRm.fastq.gz",
    log:
        "assemblies/{accession}/spades.log",
    output:
        "assemblies/{accession}_contigs.fasta",
    message:
        "Spades de novo assembly for {wildcards.accession}."
    conda:
        "../envs/spades.yaml"
    params:
        outdir="assemblies/{accession}",
    threads: 1
    shell:
        "( spades.py -1 {input.fastq_r1} -2 {input.fastq_r2} -o {params.outdir} --threads 1 --careful && "
        "mv {params.outdir}/contigs.fasta {output} ) 2> {log}"


def get_assemblies(_):
    """Get the paths to the fastqc summary files"""

    samples = pd.read_csv(
        SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species"], usecols=["Sample_Acc"]
    )

    inputs = []

    for key, sam in samples.iterrows():
        assembly = f"assemblies/{sam['Sample_Acc']}/contigs.fasta"
        inputs.append(assembly)

    return inputs


rule get_assemblies:
    input:
        get_assemblies,
    output:
        "assemblies.done",
    message:
        "Spades has finished performing the de novo assemblies."
    script:
        "touch {output}"
