#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_out_genome, get_ref_idx, get_correct_metadata


rule root_treefile:
    input:
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
    output:
        temp("trees_{pathogen}/{pathogen}_{population}_{cluster}_iq_{pruned}.nwk"),
    message:
        "Rooting tree for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/rgithub.yaml"
    wildcard_constraints:
        pruned="(pruned|rooted)",
    params:
        out=lambda wildcards: get_out_genome(wildcards),
    shell:
        "Rscript scripts/root_treefile.R {input} {params.out} {output} {wildcards.pruned}"


rule get_chromosome_ref:
    input:
        ref_index=get_ref_idx,
    output:
        temp("muttui_{pathogen}/{pathogen}_{cluster}_chr.fasta"),
        temp("muttui_{pathogen}/{pathogen}_{cluster}_pos.tsv"),
    message:
        "Extracting the chromosome for from reference genome for {wildcards.pathogen} cluster {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/get_chromosome_ref.py"


rule run_muttui_def:
    input:
        aln="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec.fasta",
        tree="muttui_{pathogen}/{pathogen}_{population}_{cluster}_iq_r.nwk",
        ref="muttui_{pathogen}/{pathogen}_{cluster}_chr.fasta",
        pos="muttui_{pathogen}/{pathogen}_{cluster}_pos.tsv",
    log:
        "muttui_{pathogen}/{pathogen}_{population}_{cluster}/muttui.log",
    output:
        "muttui_{pathogen}/{pathogen}_{population}_{cluster}/all_included_mutations.csv",
    message:
        "Running MutTui to get the mutational spectrum for {wildcards.pathogen} "
        "{wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    params:
        basename="muttui_{pathogen}/{pathogen}_{population}_{cluster}",
    shell:
        "(python ~/bin/MutTui-main/muttui/muttui.py  -a {input.aln} -t {input.tree} -r {input.ref} "
        "-c {input.pos} -o {params.basename} --exclude_root_branches) &> {log}"


rule cluster_meta:
    input:
        get_correct_metadata,
    output:
        temp("muttui_{pathogen}/{pathogen}_{population}_{cluster}_tip_labels.csv"),
    message:
        "Getting the tip labels for samples in {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/cluster_meta.py"


rule run_muttui_label:
    input:
        aln="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec.fasta",
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq_pruned.nwk",
        labels="muttui_{pathogen}/{pathogen}_{population}_{cluster}_tip_labels.csv",
    log:
        "muttui_{pathogen}/{pathogen}_{population}_{cluster}_label_tr/muttui.log",
    output:
        "muttui_{pathogen}/{pathogen}_{population}_{cluster}_label_tr/all_included_mutations.csv",
    message:
        "Running MutTui to get the mutational spectrum for each host category for "
        "{wildcards.pathogen} {wildcards.population} {wildcards.cluster} (with MutTui doing the ACE)."
    conda:
        "../envs/biopython.yaml"
    params:
        basename="muttui_{pathogen}/{pathogen}_{population}_{cluster}_label_tr",
    wildcard_constraints:
        population="(cluster|population)",
    shell:
        "(MutTui run --alignment {input.aln} --tree {input.tree}  --labels {input.labels} "
        "-o {params.basename} --exclude_root_branches --all_sites) &> {log}"
