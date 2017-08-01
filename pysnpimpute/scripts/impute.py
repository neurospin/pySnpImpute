# -*- coding: utf-8 -*-

"""
Script allowing to run the imputation of a chromosome chunk in command line.
"""

import os
import argparse
import pprint

import pysnpimpute

from pysnpimpute.utils import is_file, create_logger
from pysnpimpute.imputation import impute


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("-c", "--chromosome", required=True,
                        choices=pysnpimpute.ORDERED_CHROMOSOMES,
                        help="Name of chromosome or X chromosome region.")

    parser.add_argument("-f", "--from-bp", type=int, required=True,
                        metavar="<position>")

    parser.add_argument("-t", "--to-bp", type=int, required=True,
                        metavar="<position>")

    parser.add_argument("-i", "--hap", type=is_file, required=True,
                        metavar="<path>",
                        help="Haplotypes to impute in Impute2 format.")

    parser.add_argument("-s", "--sample", type=is_file, required=True,
                        metavar="<path>",
                        help="Samples to impute in Impute2 format.")

    parser.add_argument("-p", "--ref-hap", type=is_file, required=True,
                        metavar="<path>",
                        help="Reference haplotypes in Impute2 format.")

    parser.add_argument("-l", "--ref-legend", type=is_file, required=True,
                        metavar="<path>",
                        help="Reference legend in Impute2 format.")

    parser.add_argument("-r", "--recombination-map", type=is_file,
                        required=True, metavar="<path>",
                        help="Recombination map file for the chromosome.")

    parser.add_argument("-o", "--outdir", required=True, metavar="<path>")

    # Optional arguments
    qhelp = ("Path to the list of variants to impute. By default imputation "
             "is done for all variants of the reference panel.")
    parser.add_argument("-q", "--to-impute", metavar="<path>", help=qhelp)

    parser.add_argument("-z", "--Ne", type=int, default=20000, metavar="<int>",
                        help="Impute2 '-Ne' parameter.")

    uhelp = ("Impute2 '-buffer' parameter. Length of buffer region in Kb (NOT "
             "IN BASEPAIR) to include on each side of the analysis interval.")
    parser.add_argument("-u", "--buffer-kb", type=int, default=250,
                        metavar="<int>", help=uhelp)

    ghelp = ("By default Impute2 does not allow imputation on regions of size "
             "> 7 Mb. To force imputation on bigger regions, set this option."
             "Note that the buffers, on each side of the imputation region, "
             "have be taken into account.")
    parser.add_argument("-g", "--allow-large-regions", action="store_true",
                        help=ghelp)

    ehelp = "Path to the Impute2 executable or alias if it's in $PATH."
    parser.add_argument("-e", "--impute2-exe", default="impute2",
                        metavar="<exe>", help=ehelp)

    parser.add_argument("-Z", "--no-verbose", action="store_true",
                        help="Non-verbose mode.")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)

    return kwargs


###############################################################################
# Run script

# Get command line arguments
kwargs = get_cmd_line_args()

# Create <outdir> if it does not exist
outdir = kwargs["outdir"]
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# Print arguments
pprint.pprint(kwargs)

# If verbose, create a logger
verbose = not kwargs.pop("no_verbose")
if verbose:
    chrom = kwargs["chromosome"]
    from_bp = kwargs["from_bp"]
    to_bp = kwargs["to_bp"]
    prefix = "imputation_chr%s_%i_%i" % (chrom, from_bp, to_bp)
    logger = create_logger(prefix=prefix, log_dir=outdir)
else:
    logger = None
kwargs["logger"] = logger

# Run preprocessing
try:
    impute(**kwargs)
except Exception as e:
    if logger is not None:
        logger.error(e.message)
