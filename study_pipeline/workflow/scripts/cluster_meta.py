#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"


from utilities import population_host_metadata, get_out_genome


def cluster_meta(metadata_file, output_file, wild):
    """Function that slices out the tip labels for muttui for a given baps cluster"""

    # read the df
    meta_df = population_host_metadata(metadata_file)

    # add outgroup
    outgroup = get_out_genome(wild)

    # cluster df
    cluster_df = meta_df[
        (meta_df["cluster"] == int(wild.cluster)) | (meta_df["sample"] == outgroup)
    ]

    # write to file
    cluster_df[["sample", "Trait"]].to_csv(
        output_file, sep=",", index=False, header=False
    )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    cluster_meta(
        metadata_file=snakemake.input[0],
        output_file=snakemake.output[0],
        wild=snakemake.wildcards,
    )
