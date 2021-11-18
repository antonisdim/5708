#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

from scripts.utilities import get_right_pathogen, get_ref_genome, genome_chromosome


rule ref_dict:
    input:
        "refs/{pathogen}/{accession}.fasta",
    output:
        "refs/{pathogen}/{accession}.dict",
    message:
        "Preparing the sequene dict for {wildcards.accession}."
    conda:
        "../envs/gatk.yaml"
    shell:
        "picard CreateSequenceDictionary --REFERENCE {input} --OUTPUT {output}"


rule filter_chromosome:
    input:
        bam="bt2_alignments_{pathogen}/{sample}_ref_{accession}.bam",
        fai="refs/{pathogen}/{accession}.fasta.fai",
    output:
        temp("gatk_{pathogen}/{sample}_ref_{accession}_chr.bam"),
    message:
        "Filtering sample {wildcards.sample} keeping only the chromosome."
    conda:
        "../envs/samtools.yaml"
    params:
        chr=lambda wildcards, input: genome_chromosome(input.fai),
    shell:
        "samtools view -b -o {output} {input.bam} {params.chr}"


rule read_group:
    input:
        "gatk_{pathogen}/{sample}_ref_{accession}_chr.bam",
    log:
        "gatk_{pathogen}/{sample}_ref_{accession}_chr_RG.log",
    output:
        bam=temp("gatk_{pathogen}/{sample}_ref_{accession}_chr_RG.bam"),
        bai=temp("gatk_{pathogen}/{sample}_ref_{accession}_chr_RG.bam.bai"),
    message:
        "Preparing a bam file with readgroups for {wildcards.pathogen} sample {wildcards.sample}."
    conda:
        "../envs/gatk.yaml"
    shell:
        "(picard AddOrReplaceReadGroups --INPUT {input} --OUTPUT {output.bam} "
        "--RGLB lib1 --RGPL Illumina --RGPU Sanger --RGSM {wildcards.sample} "
        "--VALIDATION_STRINGENCY LENIENT && samtools index {output.bam}) 2> {log}"


rule haplotype_caller:
    input:
        bam_file="gatk_{pathogen}/{sample}_ref_{accession}_chr_RG.bam",
        bai_idx="gatk_{pathogen}/{sample}_ref_{accession}_chr_RG.bam.bai",
        ref="refs/{pathogen}/{accession}.fasta",
        ref_dict="refs/{pathogen}/{accession}.dict",
    log:
        "gatk_{pathogen}/{sample}_ref_{accession}.log",
    output:
        gvcf="gatk_{pathogen}/{sample}_ref_{accession}.g.vcf.gz",
        tbi="gatk_{pathogen}/{sample}_ref_{accession}.g.vcf.gz.tbi",
    message:
        "Estimating the GVCF file for {wildcards.pathogen} sample {wildcards.sample}."
    conda:
        "../envs/gatk.yaml"
    threads: 1
    shell:
        "(gatk HaplotypeCaller --input {input.bam_file} --output {output.gvcf} "
        "--reference {input.ref} --min-base-quality-score 30 "
        "--output-mode EMIT_ALL_ACTIVE_SITES --sample-ploidy 1 -ERC GVCF "
        "--do-not-run-physical-phasing true) 2> {log}"


def get_cluster_gvcfs(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(f"gatk_{wildcards.pathogen}/{sample}_ref_{ref}.g.vcf.gz")

    return input_paths


def get_ref_fasta(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{wildcards.pathogen}/{ref}.fasta"


rule combine_gvcfs:
    input:
        ref=get_ref_fasta,
        gvcfs=get_cluster_gvcfs,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_cohort.log",
    output:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_cohort.g.vcf.gz",
    message:
        "Combining GVCF files for {wildcards.pathogen} "
    conda:
        "../envs/gatk.yaml"
    shell:
        "(gatk CombineGVCFs --reference {input.ref} --variant {input.gvcfs} --output {output}) 2> {log}"
