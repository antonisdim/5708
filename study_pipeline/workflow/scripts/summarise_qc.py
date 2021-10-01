#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd


def summarise_qc(input_list, output):
    """Summarise all the fastqc outputs in one file"""

    # list to store the metrics per read file
    qc_list = []

    # read each summary file
    for infile in input_list:
        df = pd.read_csv(infile, sep="\t", names=["Status", "Metric", "Sample"])
        res = df.Status.tolist()
        res.append(df.Sample.unique()[0])
        qc_list.append(res)

    # construct the qc df and then write to disk
    qc_df = pd.DataFrame(
        qc_list,
        columns=[
            "Basic Statistics",
            "Per base sequence quality",
            "Per tile sequence quality",
            "Per sequence quality scores",
            "Per base sequence content",
            "Per sequence GC content",
            "Per base N content",
            "Sequence Length Distribution",
            "Sequence Duplication Levels",
            "Overrepresented sequences",
            "Adapter Content",
            "Sample",
        ],
    )
    qc_df.to_csv(output, sep="\t", index=False, header=True)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    summarise_qc(
        input_list=snakemake.input,
        output=snakemake.output[0],
    )
