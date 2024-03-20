#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import (
    get_out_genome,
    get_ref_idx,
    genome_chromosome,
    get_aln_file,
    get_correct_metadata,
)


rule run_iq_gtr_gamma:
    input:
        get_aln_file,
    log:
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.log",
    output:
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.bionj",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.ckp.gz",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.contree",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.iqtree",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.mldist",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.parstree",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.splits.nex",
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
    message:
        "Running iqtree for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    threads: workflow.cores * 0.25
    conda:
        "../envs/iq.yaml"
    params:
        prefix="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq",
    resources:
        mem_gb=138,
    shell:
        'cluster_input=$"{wildcards.cluster}";'
        'if [[ "$cluster_input" == "1000" ]]; then'
        "   iqtree2 -s {input} -m GTR+F+G -T {threads} -t PARS -B 1000 --prefix {params.prefix} --mem {resources.mem_gb}G 2> {log};"
        "else"
        "   (const_sites=$(snp-sites -C msa_{wildcards.pathogen}/{wildcards.pathogen}_{wildcards.population}_{wildcards.cluster}_chr_aln_nrec.fasta) && "
        "   iqtree2 -s {input} -m GTR+F+G -fconst $const_sites -T {threads} "
        "   -t PARS -B 1000 --prefix {params.prefix} "
        "   --mem {resources.mem_gb}G) 2> {log};"
        "fi"


rule run_bactdating:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq.treefile",
        pop_meta=get_correct_metadata,
    log:
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_bactdating.log",
    output:
        "trees_{pathogen}/{pathogen}_{population}_{cluster}_bactdating.tree",
    message:
        "Running bactdating on the iqtree output for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/rgithub.yaml"
    params:
        out=lambda wildcards: get_out_genome(wildcards),
        genome_size=lambda wildcards: genome_chromosome(
            get_ref_idx(wildcards), size=True
        ),
    shell:
        "(Rscript scripts/run_bactdating.R {input.tree} {params.out} {params.genome_size} "
        "{input.pop_meta} {output}) &> {log}"


rule remove_outgroup:
    input:
        "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec.fasta",
    output:
        temp(
            "msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_no_out.fasta"
        ),
    message:
        "Removing outgroup from chromosome alignment for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/remove_outgroup.py"


rule run_treetime:
    input:
        tree="trees_{pathogen}/{pathogen}_{population}_{cluster}_iq_pruned.nwk",
        alignment="msa_{pathogen}/{pathogen}_{population}_{cluster}_chr_aln_nrec_no_out.fasta",
        meta=get_correct_metadata,
    log:
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}_treetime.log",
    output:
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/trace_run.log",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/sequence_evolution_model.txt",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/molecular_clock.txt",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/skyline.tsv",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/skyline.pdf",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/timetree.pdf",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/substitution_rates.tsv",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/root_to_tip_regression.pdf",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/ancestral_sequences.fasta",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/branch_mutations.txt",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/timetree.nexus",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/divergence_tree.nexus",
        "timed_trees_{pathogen}/{pathogen}_{population}_{cluster}/dates.tsv",
    message:
        "Running treetime for {wildcards.pathogen} {wildcards.population} {wildcards.cluster}."
    conda:
        "../envs/muttui.yaml"
    resources:
        timed_tree_inst=1,
    params:
        basename="timed_trees_{pathogen}/{pathogen}_{population}_{cluster}",
        clock_stdev=0.00003,
        sky_dim=25,  # it used to be 5 now 25 for more resolution
    shell:
        "(treetime --tree {input.tree} "
        "--dates {input.meta} "
        "--name-column sample "
        "--date-column Collection_Year "
        "--aln {input.alignment} "
        "--outdir {params.basename} "
        "--report-ambiguous "
        "--confidence "
        "--clock-std-dev {params.clock_stdev} "
        "--time-marginal only-final "
        "--coalescent skyline "
        "--n-skyline {params.sky_dim} "
        "--relax 1.0 0) 2> {log}"
