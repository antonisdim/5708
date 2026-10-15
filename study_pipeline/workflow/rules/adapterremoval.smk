#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


import pandas as pd


rule adapterremoval_single_end:
    input:
        fastq="raw_reads/{accession}.fastq.gz",
    log:
        "adRm/{accession}_adRm.log",
    output:
        "adRm/{accession}_adRm.fastq.gz",
    message:
        "Trimming sequencing adapters from file {input.fastq}."
    conda:
        "../envs/adapterremoval.yaml"
    threads: 4
    params:
        basename="adRm/{accession}",
    shell:
        "(AdapterRemoval"
        "   --file1 {input}"
        "   --basename {params.basename} "
        "   --gzip "
        "   --minlength 25 "
        "   --threads {threads}"
        "   --trimns && "
        " mv {params.basename}.truncated.gz {output}"
        ") 2> {log}"


rule adapterremoval_paired_end:
    input:
        fastq_r1="raw_reads/{accession}_1.fastq.gz",
        fastq_r2="raw_reads/{accession}_2.fastq.gz",
    log:
        "adRm/{accession}_adRm.log",
    output:
        fastq_r1="adRm/{accession}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{accession}_R2_adRm.fastq.gz",
    message:
        "Trimming sequencing adapters from files {input.fastq_r1} and {input.fastq_r2}."
    conda:
        "../envs/adapterremoval.yaml"
    threads: 4
    params:
        basename="adRm/{accession}",
    shell:
        "(AdapterRemoval"
        "   --file1 {input.fastq_r1}"
        "   --file2 {input.fastq_r2} "
        "   --basename {params.basename}"
        "   --gzip "
        "   --minlength 25 "
        "   --threads {threads}"
        "   --trimns &&"
        " mv {params.basename}.pair1.truncated.gz {output.fastq_r1} && "
        " mv {params.basename}.pair2.truncated.gz {output.fastq_r2} "
        ") 2> {log}"


rule adapterremoval_collapsed:
    input:
        fastq_r1="raw_reads/{accession}_1.fastq.gz",
        fastq_r2="raw_reads/{accession}_2.fastq.gz",
    log:
        "adRm/{accession}_adRm.log",
    output:
        "adRm/{accession}_collapsed_adRm.fastq.gz",
    message:
        "Trimming sequencing adapters and collapsing reads from files {input.fastq_r1} and {input.fastq_r2}."
    conda:
        "../envs/adapterremoval.yaml"
    threads: 4
    params:
        basename="adRm/{accession}",
    shell:
        "(AdapterRemoval"
        "   --file1 {input.fastq_r1} "
        "   --file2 {input.fastq_r2} "
        "   --basename {params.basename} "
        "   --gzip "
        "   --collapse-deterministic "
        "   --minlength 25 "
        "   --threads {threads}"
        "   --trimns "
        "   --qualitymax 93 && "
        " cat {params.basename}.collapsed.gz {params.basename}.collapsed.truncated.gz 1> {output} && "
        " unlink {params.basename}.collapsed.gz && "
        " unlink {params.basename}.collapsed.truncated.gz "
        ") 2> {log}"


def get_inputs_for_count_fastq_len(wildcards):

    if (config["sample_file"] is not None) and (os.path.isfile(config["sample_file"])):
        target_samples = pd.read_csv(
            config["sample_file"],
            sep="\t",
            names=["sample_name", "read_mode", "read_file_count"],
        )
        if (
            target_samples.loc[
                target_samples["sample_name"] == wildcards.sample, "read_mode"
            ].item()
            == "PE"
        ):
            return f"adRm/{wildcards.sample}_R1_adRm.fastq.gz"
        else:
            return f"adRm/{wildcards.sample}_adRm.fastq.gz"
    else:
        if config["read_mode"] == "PE":
            return f"adRm/{wildcards.sample}_R1_adRm.fastq.gz"
        elif config["read_mode"] == "SE":
            return f"adRm/{wildcards.sample}_adRm.fastq.gz"
        else:
            return f"adRm/{wildcards.sample}_collapsed_adRm.fastq.gz"


rule count_fastq_length:
    input:
        fastq=get_inputs_for_count_fastq_len,
    output:
        "adRm/{sample}_{read_mode}.size",
    message:
        "Counting the number of reads in sample {wildcards.sample}."
    shell:
        "seqtk seq -A {input.fastq} | grep -v '^>' | wc -l 1> {output}"
