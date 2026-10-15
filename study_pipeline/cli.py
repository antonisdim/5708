#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import datetime
import os
import shutil
import sys
from multiprocessing import cpu_count
from workflow.scripts.utilities import read_sample_list, BoolType

import snakemake
import yaml
import pandas as pd
from psutil import virtual_memory

# noinspection PyUnresolvedReferences
from snakemake.api import SnakemakeApi

# noinspection PyUnresolvedReferences
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    OutputSettings,
    ResourceSettings,
    StorageSettings,
    WorkflowSettings,
    Batch,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

MAX_CPU = cpu_count()
MAX_MEM_MB = int(virtual_memory().total / 1024**2)
CONFIG_DEFAULT = "config/config.yaml"

CODE_DIR = os.path.abspath(os.path.dirname(__file__))
SNAKE_DIR = ".snakemake"

COMMANDS = ["config", "adapters", "assemblies", "alignments"]

# number to times to retry a rule that failed the first time
RESTART_TIMES = 2


class Pathomancer(object):
    """
    Command-line interface for running `pathomancer`
    """

    def __init__(self):
        parser = argparse.ArgumentParser(
            usage="""pathomancer <command> [<args>]

The pathomancer commands are:
   config       General configuration
   adapters     Remove adapters from raw illumina reads
   assemblies   Assemble samples
   alignemnts   Align reads against reference genomes


""",
        )
        parser.add_argument("command", choices=COMMANDS, help="Command to run")

        # get the CLI arguments
        args = self._parse_args(parser, level=1)

        # load the config files
        self._load_config()

        try:
            # call the class method with the name given by `command`
            reval = getattr(self, args.command)()
            exit(reval)

        except ValueError:  # todo remove
            print(ValueError)

    @staticmethod
    def _parse_args(parser, level=1):
        """
        Handle the nested-subcommand CLI interface.
        """

        # if we don't have enough arguments, then print the help
        if len(sys.argv) <= level:
            parser.print_help()
            parser.exit()

        if {"-v", "--version"}.intersection(sys.argv):
            #            print(__version__) todo bring the versioning back when I package it
            print("1.0.0")
            parser.exit()

        # slice the CLI arguments
        argv = sys.argv[1:2] if level == 1 else sys.argv[level:]

        return parser.parse_args(argv)

    def _load_config(self):  # todo think about what I want the config to include
        """
        Load the default and user space config files.
        """

        # load the default config
        with open(CONFIG_DEFAULT) as fin:
            self.config_default = yaml.safe_load(fin)

            # resolve relative paths
            self.config_default["cache"] = os.path.abspath(
                os.path.expanduser(self.config_default["cache"])
            )

        # todo not sure I want users to configure that - not sure who the users are in this instance
        # try:
        #     # load the user config
        #     with open(CONFIG_USER) as fin:
        #         self.config_user = yaml.safe_load(fin)
        #
        # except FileNotFoundError:
        #     # make the config folder (if necessary)
        #     os.makedirs(os.path.dirname(CONFIG_USER), exist_ok=True)
        #
        #     self.config_user = dict()

        self.config_user = dict()

        # merge the config dictionaries (giving precedence to the config_user)
        self.config_merged = {**self.config_default, **self.config_user}

    def adapters(self):
        """
        Build a database of target species
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="pathomancer adapter trimming",
            description="Trim sequencing adapters from raw short reads",
            add_help=False,
        )

        choice = parser.add_argument_group("Required choice")

        choice.add_argument(
            "--sample-name",
            help="Name of the sample",
            metavar="<string>",
        )

        choice.add_argument(
            "--read-mode",
            help="Paired or single end reads",
            metavar="<PE/SE>",
            # dest="db_output",
            # type=WritablePathType(),
        )

        choice.add_argument(
            "--sample-file",
            help="Tab delimited file containing sample information",
            metavar="<path>",
        )

        # add the common arguments
        self._common_arguments(parser)

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # add all command line options to the merged config
        config = {**self.config_merged, **vars(args)}

        target_list = list()

        if (config["sample_file"] is not None) and (
            os.path.isfile(config["sample_file"])
        ):
            # todo add a validation check for input file
            target_samples = read_sample_list(config["sample_file"])
            target_list = [
                f"adRm/{s}_{m}.size"
                for s, m in zip(
                    target_samples["Sample_Acc"], target_samples["Read_mode"]
                )
            ]

        elif (config["sample_name"] is not None) and (config["read_mode"] is not None):
            target_list.append(
                f"adRm/{config['sample_name']}_{config['read_mode']}.size"
            )

        else:
            raise ValueError("Invalid flag combo")

        # run the `pathomancer` workflow
        exit_code = self._run_snakemake("adapters", args, config, target_list)

        return exit_code

    def alignments(self):
        """
        Build a database of target species
        """
        # noinspection PyTypeChecker
        parser = argparse.ArgumentParser(
            prog="pathomancer read alignment",
            description="Align reads against reference genomes",
            add_help=False,
        )

        choice = parser.add_argument_group("Required choice")

        choice.add_argument(
            "--sample-name",
            help="Name of the sample",
            metavar="<string>",
        )

        choice.add_argument(
            "--read-mode",
            help="Paired or single end reads",
            metavar="<PE/SE>",
            # dest="db_output",
            # type=WritablePathType(),
        )

        choice.add_argument(
            "--species",
            help="Short code for species",
            metavar="<species>",
            # dest="db_output",
            # type=WritablePathType(),
        )

        choice.add_argument(
            "--lineage",
            help="Focus on specific lineage (lineage 1000 if the user wants to focus on species level",
            metavar="<lineage number>",
            # dest="db_output",
            # type=WritablePathType(),
        )

        choice.add_argument(
            "--sample-file",
            help="Tab delimited file containing sample information",
            metavar="<path>",
        )

        choice.add_argument(
            "--reference-genome",
            help="String OR Tab delimited file containing reference genome information",
            metavar="<path>",
        )

        optional = parser.add_argument_group("Optional arguments")

        optional.add_argument(
            "--generate-fasta",
            help=f"Generate consensus fasta pseudo-genomes (pileup with majority allele)",
            type=BoolType(),
            metavar="<bool>",
            default=argparse.SUPPRESS,
        )

        # add the common arguments
        self._common_arguments(parser)

        # get the CLI arguments
        args = self._parse_args(parser, level=2)

        # add all command line options to the merged config
        config = {**self.config_merged, **vars(args)}

        target_list = list()
        if (config["sample_file"] is not None) and (
            os.path.isfile(config["sample_file"])
        ):
            # todo add a validation check for input file
            target_samples = read_sample_list(config["sample_file"])
            target_list = [
                f"bt2_alignments_{p}/{s}_ref_{config['reference_genome']}_{m}.bam"
                for s, p, m in zip(
                    target_samples["Sample_Acc"],
                    target_samples["Species"],
                    target_samples["Read_mode"],
                )
            ]
            # if we want to generate consensus fasta sequences we extend the target list
            if config["generate_fasta"]:
                fasta_list = [
                    f"seqs_{p}/{s}_ref_{config['reference_genome']}_{m}.fasta"
                    for s, p, m in zip(
                        target_samples["Sample_Acc"],
                        target_samples["Species"],
                        target_samples["Read_mode"],
                    )
                ]
                target_list.extend(fasta_list)

        elif (config["sample_name"] is not None) and (config["read_mode"] is not None):
            target_list.append(
                f"bt2_alignments_{config['species']}/{config['sample_name']}_ref_{config['reference_genome']}_{config['read_mode']}.bam"
            )
            # we add a new target if we also want the fasta files
            if config["generate_fasta"]:
                target_list.append(
                    f"seqs_{config['species']}/{config['sample_name']}_ref_{config['reference_genome']}_{config['read_mode']}.fasta"
                )
        else:
            raise ValueError("Invalid flag combo")

        # run the `pathomancer` workflow
        exit_code = self._run_snakemake("alignments", args, config, target_list)

        return exit_code

    def _common_arguments(self, parser):
        """
        Add the common arguments shared by the `database`, `sample` and `analyse` commands
        """
        common = parser.add_argument_group("Common arguments")

        # add the help option manually so we can control where it is shown in the menu
        common.add_argument(
            "-h", "--help", action="help", help="Show this help message and exit"
        )
        common.add_argument(
            "-v", "--version", action="help", help="Print the version number and exit"
        )

        common.add_argument(
            "--cores",
            help="Maximum number of CPU cores to use",
            metavar="<int>",
            # type=IntRangeType(1, MAX_CPU),
            default=8,
        )

        common.add_argument(
            "--mem",
            help="Maximum memory (MB) to use",
            # type=IntRangeType(1024, MAX_MEM_MB),
            # default=MAX_MEM_MB if self.config_default["mem"] == "all" else self.config_default["mem"],
            metavar="<int>",
        )

        common.add_argument(
            "--unlock",
            help="Unlock the output directory following a crash or hard restart",
            action="store_true",
        )

        common.add_argument(
            "--debug",
            help="Enable debugging mode",
            action="store_true",
        )

        common.add_argument(
            "--snakemake",
            help="Pass additional flags to the `snakemake` scheduler.",
            metavar="'<json>'",
            # type=JsonType(),
        )

    def _run_snakemake(self, module, args, config, target_list):
        """
        Helper function for running the snakemake workflow
        """
        # print(f"HAYSTAC v {__version__}\n")
        print(f"Date: {datetime.datetime.now()}\n")

        config["module"] = module
        config["workflow_dir"] = os.path.join(CODE_DIR, "workflow")

        print("Config parameters:\n")
        print(module)
        params = config if args.debug else vars(args)

        for key, value in params.items():
            if value or args.debug:
                print(f" {key}: {value}")
        print("\n")

        if args.debug:
            print("Target files:\n")
            for target in target_list:
                print(" " + target)
            print("\n")

        # try:
        #     os.makedirs(SNAKE_DIR, exist_ok=True)
        # except PermissionError:
        #     print_error("Cannot write to the current working directory.")

        # # save the run-time config file
        # with open(CONFIG_RUNTIME, "w") as fout:
        #     yaml.safe_dump(config, fout, default_flow_style=False)

        # # get any extra snakemake params
        # smk_params = config.pop("snakemake") or {}

        # # get any rule targets for batching and convert them into a compatible object
        # exec_batch = Batch(config["batch"][0], config["batch"][1], config["batch"][2]) if config.get("batch") else None
        # if config.get("batch"):
        #     config.pop("batch")
        with SnakemakeApi(
            OutputSettings(
                verbose=args.debug,
                show_failed_logs=args.debug,
                printshellcmds=args.debug,
            ),
        ) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                storage_settings=StorageSettings(),
                resource_settings=ResourceSettings(
                    cores=int(args.cores),
                ),
                config_settings=ConfigSettings(config=config),
                deployment_settings=DeploymentSettings(
                    deployment_method={DeploymentMethod.CONDA},
                    conda_prefix=os.path.join("pathomancer_cache", "conda"),
                ),
                snakefile=os.path.join("workflow/workflow.smk"),
            )

            dag_api = workflow_api.dag(
                dag_settings=DAGSettings(
                    targets=target_list,
                    force_incomplete=True,
                ),
            )

            dag_api.execute_workflow(
                execution_settings=ExecutionSettings(
                    retries=0 if args.debug else RESTART_TIMES,
                    keep_going=(not args.debug),
                ),
            )

        # tidy up all the snakemake metadata
        # if len(os.listdir(os.path.join(SNAKE_DIR, "locks"))) == 0:
        #     shutil.rmtree(SNAKE_DIR)

        # translate "success" into shell exit code of 0
        # return 0 if success else 1


if __name__ == "__main__":
    Pathomancer()
