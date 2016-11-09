# -*- coding: utf-8 -*-

"""
Script allowing to run the preprocessing in command line.
"""

import os
import argparse
import pprint

from pysnpimpute.utils import create_logger
from pysnpimpute.preprocessing import preprocess

DOC = """
Preprocess the <bfile> dataset.

Steps:
- replace genders with imputed genders if requested (Plink --impute-sex)
- Remove duplicated variants (Plink --list-duplicate-vars and --exclude)
- Remove variants for which any of the following criteria is true:
    - minor allele frequency < <min_maf>
    - variant call rate < <min_vcallrate>
    - Hardy-Weinberg Equilibrium < <min_hwe>
- Remove samples for which any of the following criteria is true:
    - sample genotype call rate < <min_scallrate>
    - sex is missing (i.e. 0)
    - imputed sex differ from reported sex, if using reported sex
- lift the data to the desired genomic build, if requested. Required if
  the dataset to phase.impute is not already in the same build as the
  reference panel. Done using pyliftover.
"""


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser(
        description=DOC, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required arguments

    ihelp = "Path without extension to the dataset in PLINK BED format."
    parser.add_argument("-i", "--bfile", required=True, metavar="<path>",
                        help=ihelp)

    parser.add_argument("-o", "--outdir", required=True, metavar="<path>")

    # Optional arguments
    lhelp = ("To change the build version (e.g. from hg18 to hg19), set the "
             "source and target builds. If the data to phase/impute is not "
             "already in the same build as the reference panel, lifting is "
             "required.")
    parser.add_argument("-l", "--lift", type=str, nargs=2, metavar="<build>",
                        help=lhelp)

    thelp = "Replace reported sex with imputed sex (imputed with Plink)."
    parser.add_argument("-p", "--use-imputed-sex", action="store_true",
                        help=thelp)

    parser.add_argument("-a", "--min-maf", type=float, default=0.01,
                        metavar="<float 0-0.5>",
                        help="Minimum minor allele frequency.")

    khelp = "Minimum genotyping rate for a variant to be kept."
    parser.add_argument("-k", "--min-vcallrate", type=float, default=0.95,
                        metavar="<float 0-1>", help=khelp)

    mhelp = "Minimum genotyping rate for a subject/sample to be kept."
    parser.add_argument("-m", "--min-scallrate", type=float, default=0.95,
                        metavar="<float 0-1>", help=mhelp)

    parser.add_argument("-w", "--min-hwe", type=float, default=10e-6,
                        metavar="<float>",
                        help="Minimum Hardy-Weinberg equilibrium.")

    exe_help = "Path to the %s executable or alias if it's in $PATH."
    parser.add_argument("-e", "--plink-exe", default="plink",
                        metavar="<path>", help=exe_help % "Plink")

    parser.add_argument("-z", "--no-verbose", dest="verbose",
                        action="store_false")

    # Create a dict of arguments to use as **kwargs
    args = parser.parse_args()
    kwargs = vars(args)

    kwargs["from_build"], kwargs["to_build"] = kwargs.pop("lift")

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
verbose = kwargs.pop("verbose")
if verbose:
    logger = create_logger(prefix="preprocessing", log_dir=outdir)
else:
    logger = None
kwargs["logger"] = logger

# Run preprocessing
try:
    preprocess(**kwargs)
except Exception as e:
    if logger is not None:
        logger.error(e.message)
