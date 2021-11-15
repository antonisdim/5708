#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from functools import reduce


def summarise_alignments(aln, cov, fasta, summary_out):
    """Function to summarise alignment and coverage stats"""

    # read input files
    aln_stats = pd.read_csv(aln, sep="\t")
    cov_stats = pd.read_csv(cov, sep="\t")
    fasta_stats = pd.read_csv(fasta, sep="\t", names=["sample", "chromosome", "% N"])
    fasta_stats_pivot = fasta_stats.pivot(
        index="sample", columns="chromosome", values="% N"
    ).reset_index()
    fasta_stats_pivot.columns = [
        str(col) + "_n%" if str(col) != "sample" else str(col)
        for col in fasta_stats_pivot.columns
    ]

    # merge the iteratively
    dfs = [aln_stats, cov_stats, fasta_stats_pivot]
    summary = reduce(lambda left, right: pd.merge(left, right, on="sample"), dfs).round(
        2
    )

    # write output
    summary.to_csv(summary_out, sep="\t", index=False, header=True)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    summarise_alignments(
        aln=snakemake.input.aln,
        cov=snakemake.input.coverage,
        fasta=snakemake.input.fasta,
        summary_out=snakemake.output[0],
    )
