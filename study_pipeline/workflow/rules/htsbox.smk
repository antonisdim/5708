#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule consensus_fasta:
    input:
        bam_file="bt2_alignments_{pathogen}/{sample}_ref_{accession}.bam",
        ref="refs/{pathogen}/{accession}.fasta",
        faidx="refs/{pathogen}/{accession}.fasta.fai",
    log:
        "seqs_{pathogen}/{sample}_ref_{accession}.log",
    output:
        "seqs_{pathogen}/{sample}_ref_{accession}.fasta",
    message:
        "Creating the consensus fasta sequence for sample {wildcards.sample} for "
        "{wildcards.pathogen}."
    conda:
        "../envs/htsbox.yaml"
    shell:
        "(htsbox pileup -f {input.ref} -l 15 -T 3 -q 30 -Q 30 -M -s 3 {input.bam_file} 1> "
        "{output}) 2> {log}"
