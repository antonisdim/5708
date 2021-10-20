#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from scripts.utilities import get_ecoli_sb27, get_ecoli_all


def get_sb27_ecoli(_):
    """Get the paths for the SB27 ecoli only."""

    samples_list = get_ecoli_sb27()
    inputs_sb27_ecoli = []

    for sample in samples_list:
        inputs_sb27_ecoli.append(f"assemblies/{sample}_scaffolds.fasta")

    return inputs_sb27_ecoli


rule create_poppunk_ecoli_qfile_sb27:
    input:
        sample_list=get_sb27_ecoli,
    output:
        "poppunk_ecoli/SB27_ecoli_scaffold_paths.txt",
    message:
        "Creating the query file for PopPUNK for the SB27 E. coli samples only."
    shell:
        "for path in {input.sample_list}; do readlink -f $path; done 1> {output}"


def get_all_ecoli(_):
    """Get the paths for the SB27 ecoli plus context samples."""

    samples_list = get_ecoli_all()
    inputs_all_ecoli = []

    for sample in samples_list:
        inputs_all_ecoli.append(f"assemblies/{sample}_scaffolds.fasta")

    return inputs_all_ecoli


rule create_poppunk_ecoli_qfile_all:
    input:
        sample_list=get_all_ecoli,
    output:
        "poppunk_ecoli/SB27_plus_contx_ecoli_scaffold_paths.txt",
    message:
        "Creating the query file for PopPUNK for the SB27 and context E. coli samples."
    shell:
        "for path in {input.sample_list}; do readlink -f $path; done 1> {output}"


# #todo paths above are wrong change them
#
# rule run_poppunk_ecoli:
#     input:
#         "poppunk_ecoli/qfile.txt",
#     log:
#         "poppunk_ecoli/ecoli_cluster.log",
#     output:
#         "who knows"
#     message:
#         "Running clustering with PopPUNK for E. coli samples."
#     conda:
#         "../envs/poppunk.yaml"
#     shell:
