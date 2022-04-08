#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_right_pathogen, get_ref_genome, genome_chromosome


rule r_github:
    log:
        "aux_files/r_github.log",
    output:
        "aux_files/r_github.done",
    message:
        "Installing R packages that are not hosted on Conda."
    conda:
        "../envs/rgithub.yaml"
    shell:
        "(Rscript scripts/r_github.R &> {log}"
