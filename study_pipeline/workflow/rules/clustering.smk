#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

SAMPLE_TABLE = "samples.tsv"


# def get_ecoli_samples(_):
#     """Get the E. coli samples"""
#
#     samples = pd.read_csv(
#         SAMPLE_TABLE,
#         sep="\t",
#         names=["Sample_Acc", "Species"],
#         usecols=["Sample_Acc", "Species"],
#     )
#
#     ecoli_samples = samples[samples["Species"] == "Escherichia coli"]
#
#     inputs = []
#
#     for key, sam in ecoli_samples.iterrows():
#         inputs.append(f"adRm/{sam['Sample_Acc']}_R1_adRm.fastq.gz")
#
#     return inputs
# rule create_poppunk_ecoli_qfile:
#     input:
#         qc_check="qc/fastqc_summary.tsv",
#         sample_list=get_ecoli_samples,
#     output:
#         "poppunk_ecoli/qfile.txt",
#     message:
#         "Creating the query file for poppunk for E. coli."
#     shell:
#         "for i in {input.sample_list}; do basename $i _R1_adRm.fastq.gz | awk -F ',' "
#         '\'{{print  $1"\t" "adRm/"$1"_R1_adRm.fastq.gz\t" "adRm/"$1"_R2_adRm.fastq.gz"}}\'; '
#         "done 1> {output}"
#
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
