import pandas as pd

SAMPLE_TABLE = "samples.tsv"


def read_sample_list():
    """Read the user sample table"""

    samples = pd.read_csv(
        SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species", "Data_source"],
    )

    return samples


def get_ecoli_sb27():
    """Get the E. coli samples"""

    samples = read_sample_list()
    ecoli_samples_sb27 = samples[
        (samples["Species"] == "Escherichia coli") & (samples["Data_source"] == "SB27")
    ]
    inputs_sb27_ecoli = []

    for key, sam in ecoli_samples_sb27.iterrows():
        inputs_sb27_ecoli.append(sam['Sample_Acc'])

    return inputs_sb27_ecoli


def get_ecoli_all():
    """Get the E. coli samples"""

    samples = read_sample_list()
    ecoli_samples_all = samples[(samples["Species"] == "Escherichia coli")]
    inputs_all_ecoli = []

    for key, sam in ecoli_samples_all.iterrows():
        inputs_all_ecoli.append(sam['Sample_Acc'])

    return inputs_all_ecoli

