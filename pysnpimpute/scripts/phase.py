# -*- coding: utf-8 -*-

"""
Script allowing to run strand alignment and phasing of a chromosome
or X region in command line.
"""


import os
import argparse
import pprint

import pysnpimpute

from pysnpimpute.utils import is_file, create_logger
from pysnpimpute.phasing import phase


DOC = ("Extract data of <chromosome> from the dataset, correct strand "
       "alignment problems, if there are any, and run phasing. Variants "
       "that cannot be aligned to reference panel are removed.")


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser(description=DOC)

    # Required arguments

    parser.add_argument("-c", "--chromosome", required=True,
                        choices=pysnpimpute.ORDERED_CHROMOSOMES,
                        help="Name of chromosome or X chromosome region.")

    parser.add_argument("-i", "--bfile", required=True, metavar="<path>",
                        help="Path without extension to the dataset in Plink "
                             "bed/bim/fam format.")

    parser.add_argument("-p", "--ref-hap", type=is_file, required=True,
                        metavar="<path>",
                        help="Reference haplotypes in Impute2 format.")

    parser.add_argument("-l", "--ref-legend", type=is_file, required=True,
                        metavar="<path>",
                        help="Reference legend in Impute2 format.")

    parser.add_argument("-s", "--ref-sample", type=is_file, required=True,
                        metavar="<path>",
                        help="Reference sample in Impute2 format.")

    parser.add_argument("-r", "--recombination-map", type=is_file,
                        required=True, metavar="<path>",
                        help="Recombination map of the chromosome.")

    parser.add_argument("-o", "--outdir", required=True, metavar="<path>")

    # Optional arguments

    bhelp = ("Genomic build version of the <bfile> data (e.g. 'hg19'). "
             "Required if <chromosome> is a X region (the positions where to "
             "split depend of the genomic build).")
    parser.add_argument("-b", "--build", metavar="<version>",
                        choices=pysnpimpute.GENOMIC_BUILDS, help=bhelp)

    parser.add_argument("-z", "--Ne", type=int, default=20000, metavar="<int>",
                        help="Shapeit2 '--effective-size' parameter. "
                             "See also Impute2 '-Ne' parameter.")

    parser.add_argument("-n", "--nb-cpus", type=int, metavar="<int>",
                        help="Number of CPUs to use in the phasing step "
                             "(Shapeit2 --thread argument).")

    exe_help = "Path to the %s executable or alias if it's in $PATH."
    parser.add_argument("-j", "--plink-exe", default="plink",
                        metavar="<exe>", help=exe_help % "Plink")

    parser.add_argument("-y", "--shapeit2-exe", default="shapeit",
                        metavar="<exe>", help=exe_help % "Shapeit2")

    parser.add_argument("-Z", "--no-verbose", action="store_true",
                        help="Non-verbose mode.")

    # Create a dict of arguments to use as **kwargs
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

# If verbose (not quiet), create a logger
verbose = not kwargs.pop("no_verbose")
if verbose:
    chrom = kwargs["chromosome"]
    logger = create_logger(prefix="phasing_chr%s" % chrom, log_dir=outdir)
else:
    logger = None
kwargs["logger"] = logger

# Run phasing
try:
    phase(**kwargs)
except Exception as e:
    if logger is not None:
        logger.error(e.message)
