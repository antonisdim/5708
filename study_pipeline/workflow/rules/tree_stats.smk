#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_out_genome


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
