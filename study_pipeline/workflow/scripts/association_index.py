#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import pandas as pd
import random
import statistics

from Bio import Phylo as p
from collections import Counter
from decimal import *

from utilities import population_host_metadata


def calculate_clade_association_index(phylogeny, pop_meta_df, real_data=True):
    """
    This function takes a phylogeny (biopython) and calculates the association index of the discrete character
    at the phylogeny tips. Returns an association index.
    """

    # The index will be increased with each internal node that is analysed
    association_index = 0.0

    # Iterate through the internal nodes in the phylogeny
    for clade in phylogeny.get_nonterminals():

        # The trait list will be filled with the trait of each tip in the clade
        trait = []
        # The tip count will be increased to the number of tips present in the clade
        tip_number = 0

        # Iterate through the tips in each clade and get the sample clade
        for tip in clade.get_terminals():

            if real_data:
                # Get the trait for that tip from the pop_meta dataframe
                sample_trait = pop_meta_df.loc[
                    pop_meta_df["sample"] == str(tip), "Trait"
                ].item()
            else:
                # For bootstraps the name of the tip is the trait
                sample_trait = str(tip)

            trait.append(sample_trait)
            tip_number += 1

        # Count the number of occurrences of each trait
        trait_number = Counter(trait)

        # Assign the maximum frequency to the number of tips with the most frequent trait
        maximum_frequency = max(trait_number.values())

        # Calculate the association index of the clade
        clade_association_index = (
            Decimal(1.0) - (Decimal(maximum_frequency) / Decimal(tip_number))
        ) / ((Decimal(2.0) ** Decimal(tip_number)) - 1)
        association_index += float(clade_association_index)

    return association_index


def calculate_phylogeny_association_index(in_args):
    """Function to perform a permutation analysis and find out the association between the topology of a tree and
    certain trait"""

    # Import the newick phylogeny
    phylogeny = p.read(in_args["tree"], "newick")

    # Read the population metadata
    pop_meta_df = population_host_metadata(in_args["meta"])

    # Root and then prune the tree
    phylogeny.root_with_outgroup(outgroup_targets=in_args["outgroup"])
    phylogeny.prune(target=in_args["outgroup"])

    # Calculate the phylogeny association index
    phylogeny_association_index = calculate_clade_association_index(
        phylogeny, pop_meta_df
    )

    # List to be filled with the trait of each tip in the phylogeny
    phylogeny_tip_trait = []

    for tip in phylogeny.get_terminals():
        sample_trait = pop_meta_df.loc[
            pop_meta_df["sample"] == str(tip), "Trait"
        ].item()
        phylogeny_tip_trait.append(sample_trait)

    # List to be filled with the association index for each bootstrap run
    bootstrap_association_index = []

    # Iterate through the bootstrap replicates
    for bootstrap in range(int(in_args["boot"])):

        # Sample the traits randomly without replacement
        bootstrap_sample = random.sample(phylogeny_tip_trait, len(phylogeny_tip_trait))
        bootstrap_phylogeny = phylogeny

        # Iterate through the tips in the bootstrap phylogeny
        for idx, phylogeny_tip in enumerate(bootstrap_phylogeny.get_terminals()):
            # Assign the tip to the ith trait
            phylogeny_tip.name = bootstrap_sample[idx]

        bootstrap_association_index.append(
            calculate_clade_association_index(
                bootstrap_phylogeny, pop_meta_df, real_data=False
            )
        )

    # The boot counter will be increased if the bootstrap has a stronger association index than the real data
    number_bootstraps = 0

    # Iterate through the bootstrap association indices
    for bootstrap_index in bootstrap_association_index:

        # Check if the bootstrap association is smaller than the real data
        if bootstrap_index <= phylogeny_association_index:
            number_bootstraps += 1

    # Calculate the proportion of bootstraps with an association index as strong as the real data
    proportion_bootstraps = float(number_bootstraps) / float(in_args["boot"])

    print(
        f"Association Index of the phylogeny = {phylogeny_association_index} \n"
        f"Median bootstrap Association Index = {statistics.median(bootstrap_association_index)} \n"
        f"P-value on the association = {proportion_bootstraps}"
    )

    # Store the scores in a dataframe and then store them into a file
    scores_df = pd.DataFrame(
        [
            [
                phylogeny_association_index,
                statistics.median(bootstrap_association_index),
                proportion_bootstraps,
            ]
        ],
        columns=[
            "Phylogeny Association Index",
            "Median bootstrap Association Index",
            "P-value on the association",
        ],
    )
    scores_df.to_csv(in_args["outfile"], sep="\t", index=False, header=True)


if __name__ == "__main__":
    # Instantiate an argument parser amd populate the arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Script to calculate the association index for a trait of interest on a phylogenetic tree and "
        "calculates the statistical significance of that association. Based on the script "
        "https://github.com/chrisruis/tree_scripts/blob/main/association_index.py by Chris Ruis.",
    )
    parser.add_argument(
        "--tree", metavar="<newick file>", help="File path to newick phylogenetic tree"
    )
    parser.add_argument("--outgroup", metavar="<str>", help="Outgroup to root the tree")
    parser.add_argument(
        "--meta",
        metavar="<file>",
        help="File path to tab delimited file with the traits for the "
        "leaves of the tree",
    )
    parser.add_argument(
        "--boot", metavar="<int>", help="Number of bootstraps", default=1000
    )
    parser.add_argument("--outfile", metavar="<file>", help="Path to the output file")

    # Parse the arguments
    args = parser.parse_args()
    arg_dict = vars(args)

    # run the function
    calculate_phylogeny_association_index(arg_dict)
