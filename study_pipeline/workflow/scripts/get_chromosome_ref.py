#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import re

from Bio import SeqIO

import pandas as pd

from utilities import genome_chromosome


def get_chromosome_ref(faidx, out_fasta, out_pos):
    """Get the chromosome fasta seq of the reference genome"""

    bacterial_chromosome = []
    chr_acc = genome_chromosome(faidx)

    # get the path to the reference genome
    input_fasta = re.sub(f".gz.fai", "", os.path.abspath(faidx))

    # read the file and then store the chromosome
    with open(input_fasta) as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if record.id == chr_acc:
                bacterial_chromosome.append(record)

    with open(out_fasta, "w") as output_handle:
        SeqIO.write(bacterial_chromosome, output_handle, "fasta")

    # create the positions file
    pos = pd.DataFrame(
        data={
            "aln": [x + 1 for x in range(len(bacterial_chromosome[0].seq))],
            "ref": [x + 1 for x in range(len(bacterial_chromosome[0].seq))],
        }
    )
    pos.to_csv(out_pos, sep="\t", header=False, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    get_chromosome_ref(
        faidx=snakemake.input.ref_index,
        out_fasta=snakemake.output[0],
        out_pos=snakemake.output[1],
    )
