#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from utilities import get_mlst_header


def count_mlst(input_file, output_counts):
    """Function to count the number of samples per MLST"""

    # get the header
    header = get_mlst_header(input_file)

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
