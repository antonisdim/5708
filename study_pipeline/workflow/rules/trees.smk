#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_right_pathogen, get_ref_genome


def get_cluster_fasta_consensus(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(f"seqs_{wildcards.pathogen}/{sample}_ref_{ref}.fasta")

    return input_paths


def get_ref_idx(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{wildcards.pathogen}/{ref}.fasta.gz.fai"


rule get_chromosome_aln:
    input:
        fastas=get_cluster_fasta_consensus,
        ref_index=get_ref_idx,
    output:
        "trees_{pathogen}/{pathogen}_cluster_{cluster}_chromosome_aln.fasta",
    message:
        "Getting the chromosome alignment of {wildcards.pathogen} cluster {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/get_chromosome_aln.py"
