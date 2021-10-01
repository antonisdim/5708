#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

SAMPLE_TABLE = 'samples.tsv'

rule fastqc_paired_end:
    input:
        fastq_r1="adRm/{accession}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{accession}_R2_adRm.fastq.gz",
    log:
        "qc/{accession}_qc.log",
    output:
        "qc/{accession}_R1_adRm_fastqc.html",
        "qc/{accession}_R1_adRm_fastqc.zip",
        directory("qc/{accession}_R1_adRm_fastqc"),
        "qc/{accession}_R2_adRm_fastqc.html",
        "qc/{accession}_R2_adRm_fastqc.zip",
        directory("qc/{accession}_R2_adRm_fastqc"),
    message:
        "Running fastqc on {input.fastq_r1} and {input.fastq_r2}."
    conda:
        "../envs/fastqc.yaml"
    params:
        outdir="qc",
    shell:
        "fastqc {input.fastq_r1} {input.fastq_r2} -o {params.outdir} --extract 2> {log}"


def get_qc_samples(_):
    """Get the list of samples to run qc on"""

    samples = pd.read_csv(SAMPLE_TABLE, sep="\t", names=["Sample_Acc"])

    inputs = []

    for key, sam in samples.iterrows():
        qc_outfile = f"qc/{sam['Sample_Acc']}_R1_adRm_fastqc.zip"
        inputs.append(qc_outfile)

    return inputs



rule run_all_fastqc:
    input:
        get_qc_samples,
    output:
        "qc.done",
    message:
        "Adapters have been removed and qc has been run for all samples.",
    shell:
        "touch {output}"
