#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from Bio import SeqIO

from utilities import get_out_genome


def get_chromosome_aln(chr_alignment, wild, out_fasta):

    """Remove the outgroup sequence from the chromosome alignment"""

    bacterial_chromosomes = []
    out_acc = get_out_genome(wild)

    # parse the alignment
    with open(chr_alignment) as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):

            # if not outgroup then append to the list
            if record.id != out_acc:
                bacterial_chromosomes.append(record)

    # write the file
    with open(out_fasta, "w") as output_handle:
        SeqIO.write(bacterial_chromosomes, output_handle, "fasta")

    print(f"Now the alignment only has {len(bacterial_chromosomes)} sequences left.")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    get_chromosome_aln(
        chr_alignment=snakemake.input[0],
        wild=snakemake.wildcards,
        out_fasta=snakemake.output[0],
    )
