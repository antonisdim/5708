#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import csv


def count_mlst(input_file, output_counts):
    """Function to count the number of samples per MLST"""

    # get the number of columns
    with open(input_file) as fin:
        reader = csv.reader(fin, delimiter="\t")
        first_row = next(reader)
        num_cols = len(first_row)

    # start the header
    st_info = ["Sample", "Scheme", "ST"]
    loci = []

    # fill in the number of mlst loci
    if len(st_info) < num_cols:
        loci = [f"Locus_{num + 1}" for num in range(num_cols - len(st_info))]

    # full header
    header = st_info + loci

    # read the input file
    mlst_table = pd.read_csv(
        input_file, sep="\t", names=header, usecols=["Sample", "ST"]
    )

    # groupby, count, sort
    mlst_counts = (
        mlst_table.groupby("ST")
        .count()
        .reset_index()
        .sort_values(by="Sample", ascending=False)
    )

    # write to file
    mlst_counts.rename(columns={"Sample": "Sample count"}).to_csv(
        output_counts, sep="\t", header=True, index=False
    )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    count_mlst(
        input_file=snakemake.input[0],
        output_counts=snakemake.output[0],
    )
