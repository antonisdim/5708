#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

SAMPLE_TABLE = "samples.tsv"


rule fastqc_paired_end:
    input:
        fastq_r1="adRm/{accession}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{accession}_R2_adRm.fastq.gz",
    log:
        "qc/{accession}_qc.log",
    output:
        temp("qc/{accession}_R1_adRm_fastqc.html"),
        temp("qc/{accession}_R1_adRm_fastqc.zip"),
        temp("qc/{accession}_R1_adRm_fastqc/summary.txt"),
        temp(directory("qc/{accession}_R1_adRm_fastqc")),
        temp("qc/{accession}_R2_adRm_fastqc.html"),
        temp("qc/{accession}_R2_adRm_fastqc.zip"),
        temp("qc/{accession}_R2_adRm_fastqc/summary.txt"),
        temp(directory("qc/{accession}_R2_adRm_fastqc")),
    message:
        "Running fastqc on {input.fastq_r1} and {input.fastq_r2}."
    conda:
        "../envs/fastqc.yaml"
    params:
        outdir="qc",
    shell:
        "fastqc {input.fastq_r1} {input.fastq_r2} -o {params.outdir} --extract 2> {log}"


def get_qc_summaries(_):
    """Get the paths to the fastqc summary files"""

    samples = pd.read_csv(
        SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species"], usecols=["Sample_Acc"]
    )

    inputs = []

    for key, sam in samples.iterrows():
        qc_outfile_r1 = f"qc/{sam['Sample_Acc']}_R1_adRm_fastqc/summary.txt"
        qc_outfile_r2 = f"qc/{sam['Sample_Acc']}_R2_adRm_fastqc/summary.txt"
        inputs.append(qc_outfile_r1)
        inputs.append(qc_outfile_r2)

    return inputs


rule summarise_qc:
    input:
        get_qc_summaries,
    output:
        "qc/fastqc_summary.tsv",
    message:
        "Adapters have been removed and qc has been run for all samples."
    script:
        "../scripts/summarise_qc.py"
