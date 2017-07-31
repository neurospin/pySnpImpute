#! /usr/bin/env python
# -*- coding: utf-8 -*
##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import time
import argparse
import pprint

import pysnpimpute

from pysnpimpute.utils import (is_file,
                               load_ref_panel_description,
                               create_phased_dataset_description)
from hopla.converter import hopla


DOC = "Use Hopla to run the phasing for all the chromosomes simultaneously."

SCRIPT = os.path.join(os.path.dirname(pysnpimpute.__file__), "scripts",
                      "phase.py")


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser(description=DOC)

    # Required arguments
    parser.add_argument("-i", "--bfile", required=True, metavar="<path>",
                        help="Path without extension to the dataset in Plink "
                             "bed/bim/fam format.")

    bhelp = ("Genomic build version of the <bfile> data (e.g. 'hg19'). "
             "Required if <chromosome> is a X region (the positions where to "
             "split depend of the genomic build).")
    parser.add_argument("-b", "--build", required=True, metavar="<version>",
                        choices=pysnpimpute.GENOMIC_BUILDS, help=bhelp)

    parser.add_argument("-r", "--ref-panel", type=is_file, required=True,
                        metavar="<path>", help=pysnpimpute.REF_PANEL_DOC)

    parser.add_argument("-o", "--outdir", required=True, metavar="<outdir>")

    parser.add_argument("-n", "--nb-processes", required=True, type=int,
                        metavar="<int>")

    parser.add_argument("-c", "--nb-cpus-per-process", required=True, type=int,
                        metavar="<int>")

    # Optional arguments
    parser.add_argument("-z", "--Ne", type=int, default=20000, metavar="<int>",
                        help="Shapeit2 '--effective-size' parameter. "
                             "See also Impute2 '-Ne' parameter.")

    exe_help = "Path to the %s executable or alias if it's in $PATH."
    parser.add_argument("-j", "--plink-exe", default="plink",
                        metavar="<path>", help=exe_help % "Plink")

    parser.add_argument("-y", "--shapeit2-exe", default="shapeit",
                        metavar="<path>", help=exe_help % "Shapeit2")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)

    return kwargs


def call_hopla(bfile, ref_panel, outdir, nb_processes, nb_cpus_per_process,
               build, Ne, plink_exe, shapeit2_exe):

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Set log file name and path
    timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")
    log_filename = "hopla_phasing_%s.log" % timestamp
    log_path = os.path.join(outdir, log_filename)

    # Load the paths to the ref_panel files
    df_ref_panel = load_ref_panel_description(ref_panel)
    chromosomes = df_ref_panel.index.tolist()

    outdirs = [os.path.join(outdir, "chr" + c) for c in chromosomes]
    iterative_args = ["c", "p", "l", "s", "r", "o"]
    status, exitcodes = hopla(SCRIPT,
                              i=bfile,
                              c=chromosomes,
                              p=df_ref_panel["ref_hap"].tolist(),
                              l=df_ref_panel["ref_legend"].tolist(),
                              s=df_ref_panel["ref_sample"].tolist(),
                              r=df_ref_panel["recombination_map"].tolist(),
                              o=outdirs,
                              n=nb_cpus_per_process,
                              b=build,
                              z=Ne,
                              j=plink_exe,
                              y=shapeit2_exe,
                              hopla_iterative_kwargs=iterative_args,
                              hopla_cpus=nb_processes,
                              hopla_logfile=log_path,
                              hopla_verbose=1)

    # Create a 'phased_dataset.txt' file in <outdir>
    create_phased_dataset_description(outdir)

    for job_name, exitcode in exitcodes.items():
        if exitcode > 0:
            pprint.pprint(status[job_name]["info"])


if __name__ == "__main__":

    kwargs = get_cmd_line_args()
    pprint.pprint(kwargs)
    call_hopla(**kwargs)
