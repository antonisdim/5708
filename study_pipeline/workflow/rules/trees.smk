#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

from scripts.utilities import (
    OUT_GENOME_TABLE,
    get_right_pathogen,
    get_ref_genome,
    get_out_genome,
)


def get_cluster_fasta_consensus(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    if os.path.isfile(OUT_GENOME_TABLE):
        outgroup = get_out_genome(wildcards)
        samples_list.append(outgroup)

    for sample in samples_list:
        input_paths.append(f"seqs_{wildcards.pathogen}/{sample}_ref_{ref}.fasta")

    return input_paths


def get_ref_idx(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{wildcards.pathogen}/{ref}.fasta.gz.fai"


# todo inherit this rule
rule get_chromosome_aln:
    input:
        fastas=get_cluster_fasta_consensus,
        ref_index=get_ref_idx,
    output:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.fasta",
    message:
        "Getting the chromosome alignment of {wildcards.pathogen} cluster {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/get_chromosome_aln.py"


# todo inherit this rule
rule remove_recombination:
    input:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.fasta",
    log:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_gubbins.log",
    output:
        base_embl="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.branch_base_reconstruction.embl",
        rec_gff="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.recombination_predictions.gff",
        rec_rmbl="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.recombination_predictions.embl",
        snp_dist="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.summary_of_snp_distribution.vcf",
        branch_stats="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.per_branch_statistics.csv",
        tree_lab_final="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.node_labelled.final_tree.tre",
        log="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.log",
        tree_final="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.final_tree.tre",
        poly_phylip="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.filtered_polymorphic_sites.phylip",
        poly_fasta="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln.filtered_polymorphic_sites.fasta",
        rec_masked_fasta="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln_nrec.fasta",
    message:
        "Running Gubbins on {wildcards.pathogen} cluster {wildcards.cluster}."
    conda:
        "../envs/gubbins.yaml"
    threads: 1
    params:
        basename="msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln",
        out=lambda wildcards: get_out_genome(wildcards),
    shell:
        "(run_gubbins.py --prefix {params.basename} --tree-builder raxml "
        "--filter-percentage 70.0 --outgroup {params.out} --threads {threads} {input} && "
        "mask_gubbins_aln.py --aln {input} --gff {output.rec_gff} --out {output.rec_masked_fasta}) 2> {log}"


rule get_chromosome_snps_aln:
    input:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln_nrec.fasta",
    output:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln_nrec_snps.fasta",
    message:
        "Getting only the polymorphic positions from the chromosome alignment of "
        "{wildcards.pathogen} cluster {wildcards.cluster}."
    conda:
        "../envs/snpsites.yaml"
    shell:
        "(snp-sites -m -o {output} {input})"


rule run_raxml_gtr_gamma:
    input:
        "msa_{pathogen}/{pathogen}_cluster_{cluster}_chr_aln_nrec_snps.fasta",
    log:
        "trees_{pathogen}/RAxML_{pathogen}_cluster_{cluster}.log",
    output:
        best_tree="trees_{pathogen}/RAxML_bestTree.{pathogen}_cluster_{cluster}",
        bipartition_labels="trees_{pathogen}/RAxML_bipartitionsBranchLabels.{pathogen}_cluster_{cluster}",
        bipartition="trees_{pathogen}/RAxML_bipartitions.{pathogen}_cluster_{cluster}",
        bootstrap="trees_{pathogen}/RAxML_bootstrap.{pathogen}_cluster_{cluster}",
        info="trees_{pathogen}/RAxML_info.{pathogen}_cluster_{cluster}",
    message:
        "Running raxml for {wildcards.pathogen} cluster {wildcards.cluster}."
    threads: workflow.cores
    conda:
        "../envs/raxml.yaml"
    params:
        basename="{pathogen}_cluster_{cluster}",
        workdir=lambda wildcards: os.path.abspath(f"trees_{wildcards.pathogen}"),
    shell:
        "(raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -T {threads} -m GTRGAMMA -k -# 100 "
        "-s {input} -n {params.basename} -w {params.workdir}) 2> {log}"
