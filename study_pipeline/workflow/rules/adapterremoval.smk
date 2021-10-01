#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


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
    threads: 2
    params:
        basename="adRm/{accession}",
    shell:
        "(AdapterRemoval"
        "   --file1 {input.fastq_r1}"
        "   --file2 {input.fastq_r2} "
        "   --basename {params.basename}"
        "   --gzip "
        "   --minlength 15 "
        "   --threads {threads}"
        "   --trimns &&"
        " mv {params.basename}.pair1.truncated.gz {output.fastq_r1} && "
        " mv {params.basename}.pair2.truncated.gz {output.fastq_r2} "
        ") 2> {log}"
