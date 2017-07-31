#! /usr/bin/env python
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
import subprocess
import shutil

from collections import defaultdict

import pandas
import pysnpimpute

from pysnpimpute.utils import (is_file,
                               create_logger,
                               load_ref_panel_description,
                               load_phased_dataset_description,
                               chromosome_segmentation)
from pysnpimpute.exceptions import AutosomeNotImputable
from hopla.converter import hopla


DOC = ("Use Hopla to scale the imputation, according to the command line "
       "arguments.")
SCRIPT = os.path.join(os.path.dirname(pysnpimpute.__file__), "scripts",
                      "impute.py")


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser(description=DOC)

    ihelp = ("Path to a table file listing the paths to the phased dataset to "
             "impute (in Impute2 format) and the related chromosome names. "
             "The paths can be absolute or relative to the table file. "
             "File structure: One row per input file, tab-separated fields: "
             "<CHROM>\t<HAP>\t<SAMPLE> "
             "Accepted chromosome names: 1 .. 22, X_PAR1, X_nonPAR and X_PAR2")
    parser.add_argument("-i", "--phased-dataset", required=True, type=is_file,
                        metavar="<path>", help=ihelp)

    parser.add_argument("-r", "--ref-panel", type=is_file, required=True,
                        metavar="<path>", help=pysnpimpute.REF_PANEL_DOC)

    bhelp = ("Genomic build version of the phased dataset and reference panel "
             "(e.g. 'hg19'). Required to segment the chromosomes and not "
             "overlap with centromeres (the positions where to split depend "
             "on the genomic build).")
    parser.add_argument("-b", "--build", required=True, metavar="<version>",
                        choices=pysnpimpute.GENOMIC_BUILDS, help=bhelp)

    parser.add_argument("-o", "--outdir", required=True, metavar="<path>")

    parser.add_argument("-n", "--nb-processes", required=True, type=int,
                        metavar="<int>")

    # Optional arguments

    parser.add_argument("-z", "--Ne", type=int, default=20000, metavar="<int>",
                        help="Impute2 'Ne' parameter.")

    chelp = ("The chromosome is imputed by chunks of size ~<chunksize_kb> Kb. "
             "By default 5000 (i.e. 5Mb).")
    parser.add_argument("-k", "--chunksize-kb", type=int, default=5000,
                        metavar="<int>", help=chelp)

    Nhelp = ("The minimum number of known variants per chunk. A chunk will "
             "cover an interval of length <chunksize> Kb except if there is "
             "less than <N> variants in the interval. In that case the chunk "
             "is enlarged to ensure that there is at least N variants."
             "By default 200.")
    parser.add_argument("-N", type=int, default=200, metavar="<int>",
                        help=Nhelp)

    uhelp = ("Impute2 '-buffer' parameter. Length of buffer region in Kb (NOT "
             "IN BASEPAIR) to include on each side of the analysis interval.")
    parser.add_argument("-u", "--buffer-kb", type=int, default=250,
                        metavar="<int>", help=uhelp)

    ehelp = "Path to the Impute2 executable or alias if it's in $PATH."
    parser.add_argument("-e", "--impute2-exe", default="impute2",
                        metavar="<exe>", help=ehelp)

    parser.add_argument("-Z", "--no-verbose", dest="verbose",
                        action="store_false")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)

    return kwargs


def call_hopla(phased_dataset, ref_panel, build, outdir, nb_processes, Ne,
               chunksize_kb, N, buffer_kb, impute2_exe, verbose):

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Set up log files
    timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")
    hopla_log_name = "hopla_imputation_%s.log" % timestamp
    hopla_log_path = os.path.join(outdir, hopla_log_name)

    if verbose:
        logger = create_logger("imputation_%s" % timestamp, outdir)

    # Load the paths to the ref_panel files
    df_ref_panel = load_ref_panel_description(ref_panel)

    # Load the paths to the phased dataset to impute
    df_dataset = load_phased_dataset_description(phased_dataset)

    # Join the 2 DataFrame to map each haplotype file to impute to the
    # corresponding reference panel files
    df = pandas.concat([df_ref_panel, df_dataset], axis=1, join="inner")

    # Dict to accumulate the Hopla iterative arguments
    iter_args = defaultdict(list)

    # After imputation all the chunks of a chrom will be concatenated
    # Map <chromosome> -> list of paths of imputed chunks
    imputed_chunks_of_chrom = defaultdict(list)

    fields = ["hap", "sample", "ref_hap", "ref_legend", "ref_sample",
              "recombination_map"]
    for c, hap, samp, rhap, rleg, rsamp, rmap in df[fields].itertuples():

        # Divide the chromosome in chunks
        try:
            chunks, _, _ = chromosome_segmentation(chromosome=c,
                                                   known_hap=hap,
                                                   ref_legend=rleg,
                                                   build=build,
                                                   chunksize_kb=chunksize_kb,
                                                   N=N)
        except AutosomeNotImputable as e:
            if verbose:
                logger.warning(e.message)
            continue

        chrom_outdir = os.path.join(outdir, "chr{}".format(c))
        if not os.path.isdir(chrom_outdir):
            os.mkdir(chrom_outdir)

        for (from_bp, to_bp) in chunks:
            iter_args["c"] += [c]
            iter_args["f"] += [from_bp]
            iter_args["t"] += [to_bp]
            iter_args["i"] += [hap]
            iter_args["s"] += [samp]
            iter_args["p"] += [rhap]
            iter_args["l"] += [rleg]
            iter_args["r"] += [rmap]
            iter_args["o"] += [chrom_outdir]
#            iter_args["g"] += [(to_bp-from_bp) > 7*10**7]  # is large region
#            iter_args["q"] += [to_impute]

            # Infer path of output imputed chunk
            prefix = os.path.basename(hap).split(".gz")[0].split(".hap")[0]
            suffix = ".impute2.%i_%i" % (from_bp, to_bp)
            path_chunk = os.path.join(chrom_outdir, prefix + suffix)
            imputed_chunks_of_chrom[c] += [path_chunk]

    # Run: Hopla call
    status, exitcodes = hopla(SCRIPT,
                              z=Ne,
                              u=buffer_kb,
                              e=impute2_exe,
                              hopla_iterative_kwargs=iter_args.keys(),
                              hopla_cpus=nb_processes,
                              hopla_logfile=hopla_log_path,
                              hopla_verbose=int(verbose),
                              **iter_args)

    # Merge chunks: haplotypes and metadata files created along with imputation
    # except for _samples files that are all the same, keep one and remove the
    # others.
    for chrom in imputed_chunks_of_chrom:
        imputed_chunks = imputed_chunks_of_chrom[chrom]

        # Base path to merged chunks (remove .<start>_<end> extension)
        basepath_merged = os.path.splitext(imputed_chunks[0])[0]

        extensions = ["", "_info", "_info_by_sample", "_summary", "_warnings"]
        for ext in extensions:
            chunks_to_merge = [x + ext for x in imputed_chunks]
            path_merged = basepath_merged + ext

            with open(path_merged, "w") as f:
                subprocess.check_call(["cat"] + chunks_to_merge, stdout=f)
            subprocess.check_call(["rm"] + chunks_to_merge)

        # All _samples files are the same, keep a copy and remove the rest
        samples_files = [x + "_samples" for x in imputed_chunks]
        shutil.copyfile(samples_files[0], basepath_merged + ".samples")
        subprocess.check_call(["rm"] + samples_files)

        # Gzip the imputed merged haplotypes
        subprocess.check_call(["gzip", basepath_merged])

    for job_name, exitcode in exitcodes.items():
        if exitcode > 0:
            pprint.pprint(status[job_name]["info"])


if __name__ == "__main__":

    kwargs = get_cmd_line_args()
    pprint.pprint(kwargs)
    call_hopla(**kwargs)
