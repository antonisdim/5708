#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from functools import reduce


def summarise_tree_stats(
    fst_files, herit_files, treewas_files, assoc_idx_files, out_table
):
    """Function to summarise the different tree stats for each cluster/population"""

    # summarise the fsts

    fst_list = fst_files

    fst_df = pd.DataFrame(
        columns=[
            "Cluster",
            "Nei's Fst Rec",
            "Nei's Fst no Rec",
            "WC's Fst Rec",
            "WC's Fst no Rec",
        ]
    )

    for fst_file in fst_list:
        # read the fst file from each cluster
        fst_df_tmp = pd.read_csv(fst_file, sep="\t")

        # define temp df
        fst_cluster_df = pd.DataFrame(
            columns=[
                "Cluster",
                "Nei's Fst Rec",
                "Nei's Fst no Rec",
                "WC's Fst Rec",
                "WC's Fst no Rec",
            ]
        )

        # extract the name and the fst
        fst_cluster_df.loc[0, "Cluster"] = "_".join(fst_file.split("_")[:-1])
        fst_cluster_df.loc[0, "Nei's Fst Rec"] = fst_df_tmp.loc[
            fst_df_tmp["Rec state"] == "Recombination", "Nei's Fst"
        ].item()
        fst_cluster_df.loc[0, "Nei's Fst no Rec"] = fst_df_tmp.loc[
            fst_df_tmp["Rec state"] == "No Recombination", "Nei's Fst"
        ].item()
        fst_cluster_df.loc[0, "WC's Fst Rec"] = fst_df_tmp.loc[
            fst_df_tmp["Rec state"] == "Recombination", "WC's Fst"
        ].item()
        fst_cluster_df.loc[0, "WC's Fst no Rec"] = fst_df_tmp.loc[
            fst_df_tmp["Rec state"] == "No Recombination", "WC's Fst"
        ].item()

        # concatenate to the big dataframe
        fst_df = pd.concat([fst_df, fst_cluster_df], ignore_index=True)

    # summarise heritability

    herit_list = herit_files

    herit_df = pd.DataFrame(
        columns=["Cluster", "Heritability Rec", "Heritability no Rec"]
    )

    for herit_file in herit_list:
        # read the heritability file from each cluster
        herit_df_tmp = pd.read_csv(herit_file, sep="\t")

        # define temp df
        herit_cluster_df = pd.DataFrame(
            columns=["Cluster", "Heritability Rec", "Heritability no Rec"]
        )

        # extract the name and the H^2 values
        herit_cluster_df.loc[0, "Cluster"] = "_".join(herit_file.split("_")[:-1])
        herit_cluster_df.loc[0, "Heritability Rec"] = herit_df_tmp.loc[
            herit_df_tmp["Rec state"] == "Recombination", "Heritability (broad sense)"
        ].item()
        herit_cluster_df.loc[0, "Heritability no Rec"] = herit_df_tmp.loc[
            herit_df_tmp["Rec state"] == "No Recombination",
            "Heritability (broad sense)",
        ].item()

        # concatenate to the big dataframe
        herit_df = pd.concat([herit_df, herit_cluster_df], ignore_index=True)

    # summarise treewas

    treewas_list = treewas_files

    treewas_df = pd.DataFrame(
        columns=[
            "Cluster",
            "Terminal Num of Sig SNPs",
            "Terminal Sig Threshold",
            "Simultaneous Num of Sig SNPs",
            "Simultaneous Sig Threshold",
            "Subsequent Num of Sig SNPs",
            "Subsequent Sig Threshold",
        ]
    )

    for treewas_file in treewas_list:

        # read the treewas file from each cluster
        treewas_df_tmp = pd.read_csv(treewas_file, sep="\t")

        # define temp df
        treewas_cluster_df = pd.DataFrame(
            columns=[
                "Cluster",
                "Terminal Num of Sig SNPs",
                "Terminal Sig Threshold",
                "Simultaneous Num of Sig SNPs",
                "Simultaneous Sig Threshold",
                "Subsequent Num of Sig SNPs",
                "Subsequent Sig Threshold",
            ]
        )

        # extract the name and the treewas test values
        treewas_cluster_df.loc[0, "Cluster"] = "_".join(treewas_file.split("_")[:-1])
        treewas_cluster_df.loc[0, "Terminal Num of Sig SNPs"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Terminal", "Num of Sig SNPs"
        ].item()
        treewas_cluster_df.loc[0, "Terminal Sig Threshold"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Terminal", "Sig Threshold"
        ].item()
        treewas_cluster_df.loc[0, "Simultaneous Num of Sig SNPs"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Simultaneous", "Num of Sig SNPs"
        ].item()
        treewas_cluster_df.loc[0, "Simultaneous Sig Threshold"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Simultaneous", "Sig Threshold"
        ].item()
        treewas_cluster_df.loc[0, "Subsequent Num of Sig SNPs"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Subsequent", "Num of Sig SNPs"
        ].item()
        treewas_cluster_df.loc[0, "Subsequent Sig Threshold"] = treewas_df_tmp.loc[
            treewas_df_tmp["TreeWAS test"] == "Subsequent", "Sig Threshold"
        ].item()

        # concatenate to the big dataframe
        treewas_df = pd.concat([treewas_df, treewas_cluster_df], ignore_index=True)

    # summarise association index

    assoc_idx_list = assoc_idx_files

    assoc_idx_df = pd.DataFrame(
        columns=[
            "Cluster",
            "Phylogeny Association Index",
            "Median bootstrap Association Index",
            "P-value on the association",
        ]
    )

    for assoc_idx_file in assoc_idx_list:
        # read the fst file from each cluster
        assoc_idx_df_tmp = pd.read_csv(assoc_idx_file, sep="\t")

        # extract the name and the H^2 values
        assoc_idx_df_tmp.loc[0, "Cluster"] = "_".join(assoc_idx_file.split("_")[:-2])

        # concatenate to the big dataframe
        assoc_idx_df = pd.concat([assoc_idx_df, assoc_idx_df_tmp], ignore_index=True)

    # merge the dfs for each metric
    summary = reduce(
        lambda left, right: pd.merge(left, right, on=["Cluster"], how="outer"),
        [fst_df, herit_df, assoc_idx_df, treewas_df],
    )

    # if any fst is negative convert it to 0
    summary.set_index("Cluster", inplace=True)
    summary.clip(lower=0, inplace=True)
    summary.reset_index(inplace=True)

    # store them into a file
    summary.to_csv(out_table, sep="\t", index=False, header=True)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    summarise_tree_stats(
        fst_files=snakemake.input.fst_files,
        herit_files=snakemake.input.herit_files,
        treewas_files=snakemake.input.treewas_files,
        assoc_idx_files=snakemake.input.assoc_idx_files,
        out_table=snakemake.output[0],
    )
