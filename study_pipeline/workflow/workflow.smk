#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

from scripts.utilities import REF_GENOME_TABLE, OUT_GENOME_TABLE


include: "rules/adapterremoval.smk"
include: "rules/qc.smk"
include: "rules/clustering.smk"
include: "rules/assembly.smk"
include: "rules/mlst.smk"
include: "rules/annotation.smk"
include: "rules/pangenome.smk"
include: "rules/amr.smk"


if os.path.isfile(REF_GENOME_TABLE) or config["reference_genome"] is not None:

    include: "rules/bt2_aln.smk"
    include: "rules/htsbox.smk"
    include: "rules/seq_stats.smk"
    include: "rules/msa.smk"
    include: "rules/trees.smk"
    include: "rules/snp_calling.smk"
    include: "rules/serotyping.smk"


if os.path.isfile(OUT_GENOME_TABLE) or config["outgroup_genome"] is not None:

    include: "rules/outgroup.smk"
    include: "rules/tree_stats.smk"
    include: "rules/trait_comparisons.smk"
    include: "rules/muttui.smk"
