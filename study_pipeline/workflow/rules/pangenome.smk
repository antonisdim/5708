#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from scripts.utilities import get_right_pathogen


def get_gff_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(f"bakta_{wildcards.pathogen}/{sample}.gff3")

    return input_paths


checkpoint batch_panaroo:
    input:
        sample_list=get_gff_paths,
    output:
        "panaroo_{pathogen}/{pathogen}_batches.tsv",
    message:
        "Spltting the {wildcards.pathogen} samples into batches for panaroo."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/pangenome_batch.py"


rule convert_refseq_to_prokka_gff:
    input:
        "bakta_{pathogen}/{accession}.gff3",
    output:
        temp("bakta_{pathogen}/{accession}_prokka.gff3"),
    message:
        "Converting the bakta gff3 file for sample {wildcards.accession} into a prokka gff file."
    conda:
        "../envs/panaroo.yaml"
    shell:
        "python scripts/convert_refseq_to_prokka_gff.py --gff {input} --out {output}"


def get_gff_batch(wildcards):
    """Get the paths for the gff files."""

    batch_file = checkpoints.batch_panaroo.get(pathogen=wildcards.pathogen)
    batches = pd.read_csv(batch_file.output[0], sep="\t", names=["sample", "batch"])

    sample_batch = batches[batches["batch"] == int(wildcards.batch)]

    input_paths = []

    for idx, row in sample_batch.iterrows():
        input_paths.append(f"bakta_{wildcards.pathogen}/{row['sample']}_prokka.gff3")

    return input_paths


rule run_panaroo:
    input:
        get_gff_batch,
    log:
        "panaroo_{pathogen}/pangenome_{batch}.log",
    output:
        batch_file="panaroo_{pathogen}/pangenome_{batch}.tsv",
        outdir=directory("panaroo_{pathogen}/pangenome_{batch}"),
        pangenome_ref="panaroo_{pathogen}/pangenome_{batch}/pan_genome_reference.fa",
    message:
        "Running panaroo for pathogen {wildcards.pathogen} for batch {wildcards.batch}"
    threads: workflow.cores * 0.25
    wildcard_constraints:
        batch="\d+",
    conda:
        "../envs/panaroo.yaml"
    shell:
        '(echo {input} | awk \'BEGIN{{FS=" "; OFS="\\n"}} {{$1=$1}} 1\' > {output.batch_file} && '
        "panaroo -i {output.batch_file} -o {output.outdir} --clean-mode strict -t {threads}) 2> {log}"


def get_panaroo_graphs(wildcards):
    """Function to get the paths to the small panaroo graphs"""

    batch_file = checkpoints.batch_panaroo.get(pathogen=wildcards.pathogen)
    batches = pd.read_csv(batch_file.output[0], sep="\t", names=["sample", "batch"])
    num_of_batches = batches["batch"].max()

    return expand(
        "panaroo_{pathogen}/pangenome_{chunk_num}",
        pathogen=wildcards.pathogen,
        chunk_num=[x + 1 if num_of_batches > 1 else 1 for x in range(num_of_batches)],
    )


rule merge_panaroo_graphs:
    input:
        get_panaroo_graphs,
    log:
        "panaroo_{pathogen}/panaroo_merged.log",
    output:
        pangenome_ref="panaroo_{pathogen}/pangenome_merged/pan_genome_reference.fa",
    message:
        "Merging the panaroo graphs into a single one."
    threads: workflow.cores
    conda:
        "../envs/panaroo.yaml"
    params:
        outdir=directory("panaroo_{pathogen}/pangenome_merged"),
    shell:
        "(export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && "
        "panaroo-merge -d {input} -o {params.outdir} -t {threads}) 2> {log}"


rule get_core_msa:
    input:
        pangenome_ref="panaroo_{pathogen}/pangenome_merged/pan_genome_reference.fa",
    log:
        "panaroo_{pathogen}/panaroo_msa.log",
    output:
        msa="panaroo_{pathogen}/pangenome_merged/core_gene_alignment.aln",
    message:
        "Getting the core gene MSA for the {wildcards.pathogen} pangenome."
    threads: workflow.cores
    conda:
        "../envs/panaroo.yaml"
    params:
        outdir="panaroo_{pathogen}/pangenome_merged",
        core_thresh=0.99,
    shell:
        "(export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && "
        "panaroo-msa -o {params.outdir} -a core --aligner mafft "
        "--core_threshold {params.core_thresh} -t {threads}) 2> {log}"


rule snp_core_aln:
    input:
        msa="panaroo_{pathogen}/pangenome_merged/core_gene_alignment.aln",
    log:
        "panaroo_{pathogen}/panaroo_core_snps.log",
    output:
        snp_msa="panaroo_{pathogen}/pangenome_merged/snps_core_gene_alignment.aln",
    message:
        "Getting the SNP positions from the core gene MSA of the {wildcards.pathogen} pangenome."
    conda:
        "../envs/snpsites.yaml"
    shell:
        "(snp-sites -m -o {output.snp_msa} {input.msa} && sed -i 's/_prokka//g' {output.snp_msa}) 2> {log}"
