#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import re
import pandas as pd

from Bio import SeqIO

from utilities import get_ref_genome, genome_chromosome


def get_chromosome_aln(fasta_paths, faidx, wild, out_fasta):
    """Get the chromosome fasta seq"""

    bacterial_chromosomes = []
    chr_acc = genome_chromosome(faidx)
    ref_acc = get_ref_genome(wild)

    for fasta_file in fasta_paths:

        # get sample name
        sample_name = re.sub(f"_ref_{ref_acc}.fasta", "", os.path.basename(fasta_file))

        # parse the file and get the chromosome sequence (the largest one)
        for record in SeqIO.parse(fasta_file, "fasta"):

            if record.id == chr_acc:
                record.id = sample_name
                record.name = ""
                record.description = ""
                bacterial_chromosomes.append(record)

            print(f"Done with sample {sample_name}.")

    with open(out_fasta, "w") as output_handle:
        SeqIO.write(bacterial_chromosomes, output_handle, "fasta")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    get_chromosome_aln(
        fasta_paths=snakemake.input.fastas,
        faidx=snakemake.input.ref_index,
        wild=snakemake.wildcards,
        out_fasta=snakemake.output[0],
    )
