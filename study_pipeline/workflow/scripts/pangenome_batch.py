__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import os

import random

BATCH_SIZE = 500
SEED = 12345


def pangenome_batch(gff_paths, output):
    """Split gff files into batches of 500 files"""

    random.Random(SEED).shuffle(gff_paths)

    batch = 1
    counter = 0
    batch_list = []

    for gff in gff_paths:

        sample_name = os.path.splitext(os.path.basename(gff))[0]

        if counter >= BATCH_SIZE:
            counter = 1
            batch += 1
        else:
            counter += 1

        batch_list.append([sample_name, batch])

    batch_df = pd.DataFrame(batch_list, columns=["sample", "batch"])
    batch_df.to_csv(output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    pangenome_batch(
        gff_paths=snakemake.input.sample_list,
        output=snakemake.output[0],
    )
