#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_right_pathogen


def get_gff_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(f"prokka_{wildcards.pathogen}/{sample}.gff")

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


def get_gff_batch(wildcards):
    """Get the paths for the gff files."""

    batch_file = checkpoints.batch_panaroo.get(pathogen=wildcards.pathogen)
    batches = pd.read_csv(batch_file.output[0], sep="\t", names=["sample", "batch"])

    sample_batch = batches[batches["batch"] == int(wildcards.batch)]

    input_paths = []

    for idx, row in sample_batch.iterrows():
        input_paths.append(f"prokka_{wildcards.pathogen}/{row['sample']}.gff")

    return input_paths


rule run_panaroo:
    input:
        get_gff_batch,
    log:
        "panaroo_{pathogen}/pangenome_{batch}.log",
    output:
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
        "(panaroo -i {input} -o {output.outdir} --clean-mode strict -t {threads}) 2> {log}"


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
        "panaroo_{pathogen}/pangenome_merged/panaroo.log",
    output:
        outdir=directory("panaroo_{pathogen}/pangenome_merged"),
        pangenome_ref="panaroo_{pathogen}/pangenome_merged/pan_genome_reference.fa",
    message:
        "Merging the panaroo graphs into a single one"
    threads: workflow.cores
    conda:
        "../envs/panaroo.yaml"
    shell:
        "(panaroo-merge -d {input} -o {output.outdir} -t {threads}) 2> {log}"
