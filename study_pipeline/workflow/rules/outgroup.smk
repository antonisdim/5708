#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule split_outgroup:
    input:
        "refs/{pathogen}/{accession}.fasta.gz",
    output:
        unzip=temp("refs/{pathogen}/{accession}_temp.fasta"),
        flat=temp("refs/{pathogen}/{accession}_temp.fasta.flat"),
        gdx=temp("refs/{pathogen}/{accession}_temp.fasta.gdx"),
        kmer=temp("refs/{pathogen}/{accession}_temp.split.500mer.5overlap.fasta"),
    message:
        "Splitting the outgroup reference genome {wildcards.accession} into kmers "
        "of 500 bp for {wildcards.pathogen}."
    conda:
        "../envs/pyfasta.yaml"
    shell:
        "zcat {input} > {output.unzip} && pyfasta split -n 1 -k 500 -o 5 {output.unzip}"


rule bowtie_aln_outgroup:
    input:
        fasta="refs/{pathogen}/{sample}_temp.split.500mer.5overlap.fasta",
        db_idx="refs/{pathogen}/{accession}.1.bt2l",
    log:
        "bt2_alignments_{pathogen}/{sample}_{accession}.log",
    output:
        bam_file="bt2_alignments_{pathogen}/{sample}_ref_{accession}.bam",
        bai_file="bt2_alignments_{pathogen}/{sample}_ref_{accession}.bam.bai",
    params:
        basename="refs/{pathogen}/{accession}",
    threads: 4
    message:
        "Aligning the reads from sample {wildcards.sample} against taxon {wildcards.pathogen}."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "(bowtie2 -f --time --no-unal --no-mixed --threads {threads} -x {params.basename} -U {input.fasta} | "
        "samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file}) 2> {log}"
