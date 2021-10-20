import pandas as pd

SAMPLE_TABLE = "samples.tsv"


def read_sample_list():
    """Read the user sample table"""

    samples = pd.read_csv(
        SAMPLE_TABLE, sep="\t", names=["Sample_Acc", "Species"],
    )

    return samples