#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule coverage_counts:
    input:
        bam_file="bt2_alignments_{pathogen}/{sample}_ref_{accession}.bam",
    output:
        temp("bt2_alignments_{pathogen}/{sample}_cov_count_{accession}.txt"),
    message:
        "Counting coverage stats for sample {wildcards.sample} and taxon {wildcards.pathogen}."
    conda:
        "../envs/samtools.yaml"
    shell:
        "(samtools mpileup -a {input} | "
        " awk '$4 != \"0\" {{rows++; sum += $4}} END {{print NR, rows, sum}}' OFS='\t'"
        ") 1> {output} 2> /dev/null"


rule coverage_stats:
    input:
        "bt2_alignments_{pathogen}/{sample}_cov_count_{accession}.txt",
        "refs/{pathogen}/{accession}.fasta.gz.fai",
    output:
        temp("bt2_alignments_{pathogen}/{sample}_cov_stats_{accession}.txt"),
    message:
        "Calculating coverage statistics to assess if reads from sample {wildcards.sample} "
        "represent a random genome sample of taxon {wildcards.pathogen}."
    script:
        "../scripts/coverage_stats.py"


def get_cluster_stats(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(
            f"bt2_alignments_{wildcards.pathogen}/{sample}_cov_stats_{ref}.txt"
        )

    return input_paths


rule summarise_cluster_stats:
    input:
        get_cluster_stats,
    output:
        "bt2_alignments_{pathogen}/{pathogen}_cluster_{cluster}_coverage_stats.tsv",
    message:
        "Concatenating all the coverage stats for {wildcards.pathogen}."
    shell:
        "printf "
        '"sample\ttaxon\tref_bases_cov\ttotal_bases_cov\tcoverage\tfraction_ref_cov\tcov_evenness\n" > '
        "{output} && cat {input} >> {output}"


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


rule unknown_bases:
    input:
        "seqs_{pathogen}/{sample}_ref_{accession}.fasta",
    output:
        temp("seqs_{pathogen}/{sample}_ncount_{accession}.txt")
    message:
        "Counting the number of unknown bases in the consenus fasta sequence for sample "
        "{wildcards.sample} for {wildcards.pathogen}."
    conda:
        "../envs/seq.yaml"
    shell:
        "seqtk comp {input} | column -t | awk '{{print {wildcards.sample}\"\t\"$1\"\t\"$9/$2*100}}' > "
        "{output}"


def get_cluster_ncounts(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(
            f"seqs_{wildcards.pathogen}/{sample}_ncount_{ref}.txt"
        )

    return input_paths


rule cat_unknown_base_counts:
    input:
        get_cluster_ncounts
    output:
        "seqs_{pathogen}/{pathogen}_cluster_{cluster}_n_stats.tsv"
    message:
        "Getting the number of Ns in consensus fasta seqs across all samples "
        "for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "cat {input} > {output}"