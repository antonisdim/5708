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
    get_ref_idx,
    get_aln_file,
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


rule get_chromosome_aln:
    input:
        fastas=get_cluster_fasta_consensus,
        ref_index=get_ref_idx,
    output:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.fasta",
    message:
        "Getting the chromosome alignment of {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/get_chromosome_aln.py"


rule remove_recombination:
    input:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.fasta",
    log:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_gubbins.log",
    output:
        base_embl="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.branch_base_reconstruction.embl",
        rec_gff="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.recombination_predictions.gff",
        rec_rmbl="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.recombination_predictions.embl",
        snp_dist="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.summary_of_snp_distribution.vcf",
        branch_stats="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.per_branch_statistics.csv",
        tree_lab_final="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.node_labelled.final_tree.tre",
        log="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.log",
        tree_final="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.final_tree.tre",
        poly_phylip="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.filtered_polymorphic_sites.phylip",
        poly_fasta="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.filtered_polymorphic_sites.fasta",
        rec_masked_fasta="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec.fasta",
    message:
        "Running Gubbins on {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/gubbins.yaml"
    threads: workflow.cores * 0.5
    params:
        basename="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln",
        out=lambda wildcards: get_out_genome(wildcards),
        percent_n=95.0,
    shell:
        "(export OPENBLAS_NUM_THREADS=1 && "
        "run_gubbins.py --prefix {params.basename} --first-tree-builder fasttree --tree-builder raxml "
        "--filter-percentage {params.percent_n} --outgroup {params.out} --threads {threads} {input} && "
        "mask_gubbins_aln.py --aln {input} --gff {output.rec_gff} --out {output.rec_masked_fasta}) 2> {log}"


rule get_chromosome_snps_aln:
    input:
        get_aln_file,
    output:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_{rec}_snps.fasta",
    message:
        "Getting only the polymorphic positions from the {wildcards.rec} chromosome alignment of "
        "{wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/snpsites.yaml"
    shell:
        "(snp-sites -m -o {output} {input})"
