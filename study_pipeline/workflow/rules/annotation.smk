#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_right_pathogen


rule run_prokka:
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
    shell:
        '(if [[ "{wildcards.pathogen}" == "ecoli" ]]; then '
        "   prokka "
        "       --outdir {params.outdir} "
        "       --force -prefix {wildcards.accession} "
        "       --addgenes "
        "       --centre Sanger "
        "       --compliant "
        "       --genus Escherichia --species coli "
        "       --kingdom Bacteria "
        "       --gcode 11 "
        "       --usegenus "
        "       --cpus {threads} "
        "       --evalue 1e-6 "
        "       {input}; "
        "else "
        "   prokka "
        "       --outdir {params.outdir} "
        "       --force -prefix {wildcards.accession} "
        "       --addgenes "
        "       --centre Sanger "
        "       --compliant"
        "       --kingdom Bacteria "
        "       --gcode 11 "
        "       --cpus {threads} "
        "       --evalue 1e-6 "
        "       {input}; "
        "fi) 2> {log}"


rule bakta_db:
    log:
        "public_dbs/db_bakta.log",
    output:
        directory("public_dbs/db"),
    message:
        "Downloading the Bakta database."
    conda:
        "../envs/bakta.yaml"
    params:
        outdir="public_dbs",
    shell:
        "bakta_db download --output {params.outdir} 2> {log}"


rule run_bakta:
    input:
        scaffold="assemblies/{accession}_scaffolds.fasta",
        db="public_dbs/db",
    log:
        "bakta_{pathogen}/{accession}_screen.log",
    output:
        "bakta_{pathogen}/{accession}.embl",
        "bakta_{pathogen}/{accession}.faa",
        "bakta_{pathogen}/{accession}.ffn",
        "bakta_{pathogen}/{accession}.fna",
        "bakta_{pathogen}/{accession}.gbff",
        "bakta_{pathogen}/{accession}.gff3",
        "bakta_{pathogen}/{accession}.hypotheticals.faa",
        "bakta_{pathogen}/{accession}.hypotheticals.tsv",
        "bakta_{pathogen}/{accession}.json",
        "bakta_{pathogen}/{accession}.tsv",
        "bakta_{pathogen}/{accession}.txt",
        "bakta_{pathogen}/{accession}.log",
    message:
        "Running bakta for {wildcards.pathogen} sample with acc {wildcards.accession}."
    conda:
        "../envs/bakta.yaml"
    threads: 1
    params:
        outdir="bakta_{pathogen}",
    shell:
        "(bakta --db {input.db} --output {params.outdir} --prefix {wildcards.accession} "
        "--verbose --threads {threads}  {input.scaffold}) &> {log}"


def get_bakta_gff_paths(wildcards):
    """Get the paths for the SB27 (and context samples)."""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    for sample in samples_list:
        input_paths.append(f"bakta_{wildcards.pathogen}/{sample}.gff3")

    return input_paths


rule aggregate_annotation:
    input:
        get_bakta_gff_paths,
    output:
        "annotation_{pathogen}.done",
    message:
        "All annotation is done."
    shell:
        "touch {output}"
