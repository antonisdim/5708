#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_right_pathogen, get_ref_genome

MIN_FRAG_LEN = 0
MAX_FRAG_LEN = 1000


rule samtools_index_accession:
    input:
        "refs/{pathogen}/{accession}.fasta.gz",
    output:
        "refs/{pathogen}/{accession}.fasta.gz.fai",
    message:
        "Indexing fasta file with accession {wildcards.accession} for taxon {wildcards.pathogen}."
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule bowtie_index_accession:
    input:
        "refs/{pathogen}/{accession}.fasta.gz",
    log:
        "refs/{pathogen}/{accession}_index.log",
    output:
        "refs/{pathogen}/{accession}.1.bt2l",
        "refs/{pathogen}/{accession}.2.bt2l",
        "refs/{pathogen}/{accession}.3.bt2l",
        "refs/{pathogen}/{accession}.4.bt2l",
        "refs/{pathogen}/{accession}.rev.1.bt2l",
        "refs/{pathogen}/{accession}.rev.2.bt2l",
    message:
        "Preparing the bowtie2 index for genome {wildcards.accession} of taxon {wildcards.pathogen}."
    threads: 4
    params:
        basename="refs/{pathogen}/{accession}",
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --large-index --threads {threads} {input} {params.basename} &> {log}"


rule bowtie_align_accession_paired_end:
    input:
        fastq_r1="adRm/{sample}_R1_adRm.fastq.gz",
        fastq_r2="adRm/{sample}_R2_adRm.fastq.gz",
        db_idx="refs/{pathogen}/{accession}.1.bt2l",
    log:
        "bt2_alignments_{pathogen}/{sample}_{accession}.log",
    output:
        bam_file="bt2_alignments_{pathogen}/{sample}-{accession}.bam",
        bai_file="bt2_alignments_{pathogen}/{sample}-{accession}.bam.bai",
    params:
        basename="refs/{pathogen}/{accession}",
    threads: 4
    message:
        "Aligning the reads from sample {wildcards.sample} against taxon {wildcards.pathogen}."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 --time --no-unal --no-mixed --threads {threads} "
        "-I {MIN_FRAG_LEN} -X {MAX_FRAG_LEN} -x {params.basename} "
        "-1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file} "
        ") 2> {log}"


def get_cluster_samples(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(f"bt2_alignments_{wildcards.pathogen}/{sample}-{ref}.bam")

    return input_paths


rule aln_done:
    input:
        get_cluster_samples,
    output:
        "aln_{pathogen}_lin{cluster}.done",
    message:
        "All local alignemnts with bt2 are done for cluster {wildcards.cluster} "
        "and pathogen {wildcards.pathogen}"
    shell:
        "touch {output}"
