#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import read_sample_list


rule run_prokka_ecoli:
    input:
        "assemblies/{accession}_scaffolds.fasta",
    log:
        "prokka_{pathogen}/{accession}.log",
    output:
        "prokka_{pathogen}/{accession}.err",
        "prokka_{pathogen}/{accession}.faa",
        "prokka_{pathogen}/{accession}.ffn",
        "prokka_{pathogen}/{accession}.fna",
        "prokka_{pathogen}/{accession}.fsa",
        "prokka_{pathogen}/{accession}.gbk",
        "prokka_{pathogen}/{accession}.gff",
        "prokka_{pathogen}/{accession}.sqn",
        "prokka_{pathogen}/{accession}.tbl",
        "prokka_{pathogen}/{accession}.tsv",
        "prokka_{pathogen}/{accession}.txt",
    message:
        "Running prokka for {wildcards.pathogen} samples with acc {wildcards.accession}."
    conda:
        "../envs/panaroo.yaml"
    threads: 1
    params:
        outdir="prokka_{pathogen}",
        db="{pathogen}_poppunk",
    shell:
        "(prokka"
        "   --outdir {params.outdir}"
        "   --force -prefix {wildcards.accession}"
        "   --addgenes"
        "   --centre Sanger"
        "   --compliant"
        "   --genus Escherichia --species coli"
        "   --kingdom Bacteria"
        "   --gcode 11"
        "   --usegenus"
        "   --cpus {threads}"
        "   --evalue 1e-6"
        "   {input} ) 2> {log}"


def get_gff_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(f"prokka_{wildcards.pathogen}/{sample}.gff")

    return input_paths


rule aggregate_annotation:
    input:
        get_gff_paths,
    output:
        "annotation_{pathogen}.done",
    message:
        "All annotation is done."
    shell:
        "touch {output}"
