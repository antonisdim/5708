#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import (
    get_out_genome,
    genome_chromosome,
    get_ref_idx,
    get_clock_rate,
    TIME_INTERVAL,
)


rule pairsnp:
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


# todo change the pop meta in a wild card combo that fits both clusters and pops or combine them into a single file
rule transition_analysis:
    input:
        snps="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_pairsnp.tsv",
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_transition_analysis.log",
    output:
        all_hosts="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_total_host_links.tsv",
        boot_hosts="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_boot_host_links.tsv",
        anc_states="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_ace_summary.tsv",
        anc_counts="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_ace_tr_summary.tsv",
    message:
        "Inferring host transition links for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    params:
        out=lambda wildcards: get_out_genome(wildcards),
        clock=lambda wildcards: get_clock_rate(wildcards),
        genome_size=lambda wildcards: genome_chromosome(
            get_ref_idx(wildcards), size=True
        ),
    conda:
        "../envs/rgithub.yaml"
    shell:
        "(Rscript scripts/transition_analysis.R {input.snps} {input.tree} {params.out} {input.pop_meta} "
        "{output.all_hosts} {output.boot_hosts} {output.anc_states} {output.anc_counts} {params.clock} "
        "{params.genome_size} {TIME_INTERVAL}) &> {log}"


rule cohort_comparison:
    input:
        aln="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_snps.fasta",
        pop_meta="aux_files/{pathogen}_all_meta.tsv",
    log:
        "trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_{set}_cohort_comparison.log",
    output:
        histogram="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_{set}_cohort_hist.pdf",
        pvalues="trees_stats_{pathogen}/{pathogen}_{population}_{cluster}_{set}_cohort_pvalues.tsv",
    message:
        "Comparing whether the {wildcards.set} cohort differs from the contextual/animal one "
        "for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    params:
        outgroup=lambda wildcards: get_out_genome(wildcards),
    conda:
        "../envs/rgithub.yaml"
    shell:
        "(Rscript scripts/cohort_comparison.R {params.outgroup} {input.aln} {input.pop_meta} "
        "{output.histogram} {output.pvalues} {wildcards.set}) &> {log}"
