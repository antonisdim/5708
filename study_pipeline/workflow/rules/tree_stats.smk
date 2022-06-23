#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_out_genome, get_clusters_to_run


rule github_r:
    log:
        "aux_files/github_r.log",
    output:
        "aux_files/github_r.done",
    message:
        "Installing R packages that are not hosted on Conda."
    conda:
        "../envs/rgithub.yaml"
    shell:
        "(Rscript scripts/github_r.R {output}) &> {log}"


rule treewas:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        aln="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_snps.fasta",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_treewas.log",
    output:
        plot="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_treewas.pdf",
        table="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_treewas.tsv",
    message:
        "Running treeWAS on {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/rgithub.yaml"
    params:
        out=lambda wildcards: get_out_genome(wildcards),
    shell:
        "(Rscript scripts/treewas.R {input.tree} {params.out} {input.aln} {input.pop_meta} "
        "{output.plot} {output.table}) &> {log}"


# todo move that to the original rule in trees.smk - here only for convenience so that it doesn't trigger the re-execution of trees
rule get_chromosome_snps_nrec_aln:
    input:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln.fasta",
    output:
        temp("msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_snps.fasta"),
    message:
        "Getting only the polymorphic positions from the chromosome alignment of "
        "{wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/snpsites.yaml"
    shell:
        "(snp-sites -m -o {output} {input})"


rule fst:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        aln_rec="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_snps.fasta",
        aln_nrec="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_snps.fasta",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_fst.log",
    output:
        table="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_fst.tsv",
    message:
        "Calucluate Nei's and Weir and Cockerham's Fst for {wildcards.pathogen} {wildcards.population} "
        "{wildcards.cluster}, before and after removing recombinant regions."
    conda:
        "../envs/rgithub.yaml"
    params:
        out=lambda wildcards: get_out_genome(wildcards),
    shell:
        "(Rscript scripts/fst.R {input.tree} {params.out} {input.aln_rec} {input.aln_nrec} "
        "{input.pop_meta} {output.table}) &> {log}"


rule heritability:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        aln_rec="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_snps.fasta",
        aln_nrec="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_snps.fasta",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_heritability.log",
    output:
        table="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_heritability.tsv",
    message:
        "Calucluate the broad sense heritability (H^2) from an amova, "
        "for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}, "
        "before and after removing recombinant regions."
    conda:
        "../envs/rgithub.yaml"
    params:
        out=lambda wildcards: get_out_genome(wildcards),
    shell:
        "(Rscript scripts/heritability.R {input.tree} {params.out} {input.aln_rec} {input.aln_nrec} "
        "{input.pop_meta} {output.table}) &> {log}"


rule association_index:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_assoc_idx.log",
    output:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_assoc_idx.tsv",
    message:
        "Calucluate the phylogeny and trait (host) association index, "
        "for {wildcards.pathogen} {wildcards.population} {wildcards.cluster} "
        "for {params.boot} bootstraps, for the clonal tree only."
    conda:
        "../envs/biopython.yaml"
    params:
        out=lambda wildcards: get_out_genome(wildcards),
        boot=1000,
    shell:
        "(python scripts/association_index.py --tree {input.tree} --boot {params.boot} "
        "--meta {input.pop_meta} --outgroup {params.out} --outfile {output}) 2> {log}"


def get_cluster_treewas(wildcards):
    """Get the fst files belonging for the clusters that have been run so far"""

    input_paths = []
    cluster_list = get_clusters_to_run()

    for cluster in cluster_list:
        input_paths.append(
            f"trees_stats_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{cluster}_treewas.tsv"
        )

    return input_paths


def get_cluster_fst(wildcards):
    """Get the fst files belonging for the clusters that have been run so far"""

    input_paths = []
    cluster_list = get_clusters_to_run()

    for cluster in cluster_list:
        input_paths.append(
            f"trees_stats_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{cluster}_fst.tsv"
        )

    return input_paths


def get_cluster_heritability(wildcards):
    """Get the fst files belonging for the clusters that have been run so far"""

    input_paths = []
    cluster_list = get_clusters_to_run()

    for cluster in cluster_list:
        input_paths.append(
            f"trees_stats_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{cluster}_heritability.tsv"
        )

    return input_paths


def get_cluster_assoc_idx(wildcards):
    """Get the fst files belonging for the clusters that have been run so far"""

    input_paths = []
    cluster_list = get_clusters_to_run()

    for cluster in cluster_list:
        input_paths.append(
            f"trees_stats_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{cluster}_assoc_idx.tsv"
        )

    return input_paths


rule summarise_tree_stats:
    input:
        fst_files=get_cluster_fst,
        herit_files=get_cluster_heritability,
        treewas_files=get_cluster_treewas,
        assoc_idx_files=get_cluster_assoc_idx,
    output:
        "trees_stats_{pathogen}/{pathogen}_{population}_tree_stats.tsv",
    message:
        "Summarising the tree stats for {wildcards.pathogen} {wildcards.population}(s), using host as a trait."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/summarise_tree_stats.py"


rule  pairsnp:
    input:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_snps.fasta",
    output:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_pairsnp.tsv",
    message:
        "Running pairsnp for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/pairsnp.yaml"
    shell:
        "pairsnp -i {input} > {output}"


rule transition_analysis:
    input:
        snps="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_pairsnp.tsv",
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    output:
        all_hosts="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_total_host_links.tsv",
        boot_hosts="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_boot_host_links.tsv"
    message:
        "Inferring host transition links for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    params:
        out=lambda wildcards: get_out_genome(wildcards),
    conda:
        "../envs/rgithub.yaml"
    shell:
        "(Rscript scripts/transition_analysis.R {input.snps} {input.tree} {params.out} {input.pop_meta} "
        "{output.all_hosts} {output.boot_hosts})"

