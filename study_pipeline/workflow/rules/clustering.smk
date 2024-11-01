#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


checkpoint run_fastbaps:
    input:
        "panaroo_{pathogen}/pangenome_merged/snps_core_gene_alignment.aln",
    log:
        "clustering_{pathogen}/fastbaps.log",
    output:
        "clustering_{pathogen}/fastbaps_clusters_{pathogen}.csv",
    message:
        "Running fastbaps for the SB27 and context samples of {wildcards.pathogen}."
    conda:
        "../envs/fastbaps.yaml"
    # threads: workflow.cores
    shell:
        "(Rscript scripts/run_fastbaps.R {input} {output} 1) &> {log}"


rule fastbaps_summary:
    input:
        "clustering_{pathogen}/fastbaps_clusters_{pathogen}.csv",
    output:
        "clustering_{pathogen}/summary_fastbaps_clusters_{pathogen}.tsv",
    message:
        "Quick summary of the fastbaps cluster output for the SB27 samples of {wildcards.pathogen}."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/fastbaps_summary.py"


rule clusters_and_mlst:
    input:
        clusters_summary=(
            "clustering_{pathogen}/summary_fastbaps_clusters_{pathogen}.tsv"
        ),
        mlst_summary="mlst/mlst_summary_SB27_{pathogen}_{dataset}_contx.tsv",
        clusters="clustering_{pathogen}/fastbaps_clusters_{pathogen}.csv",
        mlst="mlst/SB27_{pathogen}_{dataset}_contx.tsv",
    output:
        "clustering_{pathogen}/summary_clusters_mlst_{pathogen}_{dataset}_contx.tsv",
    message:
        "Counting the number of STs per cluster for {wildcards.pathogen}."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/clusters_and_mlst.py"
