# -*- coding: utf-8 -*-

"""
Defines a set of functions to run the imputation using Impute2.

The module requires Impute2 to be installed.
"""

import os

import pysnpimpute

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_chromosome_name,
                               check_existence_of_paths,
                               run_cmd_and_check)


def impute(chromosome, from_bp, to_bp, hap, sample, ref_hap, ref_legend,
           recombination_map, outdir=None, to_impute=None, Ne=20000,
           buffer_kb=250, allow_large_regions=False, basename=None,
           suffix=None, impute2_exe="impute2", logger=None):
    """
    Run imputation using Impute2

    Parameters
    ----------
    chromosome: str
        Name of chromosome or X chromosome region.
        Accepted names: "1", ..., "22", "X_PAR1", "X_nonPAR" and "X_PAR2".
    from_bp, to_bp: int
        The interval in basepair position to impute.
    hap: str
        Path to the phased data to impute in Impute2 format.
    sample: str
        Path to samples to impute in Impute2 format.
    ref_hap, ref_legend: str
        Path to reference panel file in Shapeit2/Impute2 format.
    recombination_map: str
        Path to the recombination map required by Shapeit2 and Impute2.
    outdir: str
        Path to directory where to output.
    to_impute: str, default None
        Path to the list of variants to impute. By default imputation is done
        for all variants of the reference panel.
    Ne: int, default 20000
        Impute2 'Ne' paraemter.
    buffer_kb: int, default 250
        Impute2 '-buffer' parameter. Length of buffer region in kb (NOT IN
        BASEPAIR) to include on each side of the analysis interval.
    allow_large_regions: bool, default False
        By default Impute2 does not allow imputation on a region of size > 7Mb.
        To force imputation on a bigger region, set this option to True.
    basename, suffix: str, default None.
        Output path is <outdir>/<basename><suffix>.
        By default basename is <hap> filename without .hap.gz extension and
        suffix is '.impute2.<from_bp>_<to_bp>'.
    impute2_exe: str, default "minimac3"
        Path to the impute2 executable or alias if it's in $PATH.
    logger: logging object, defaut None.
        To activate logging, pass a logging object.
    """

    # Check that Impute2 is installed
    check_installation_of_required_softwares(dict(Impute2=impute2_exe))
    check_chromosome_name(chromosome)

    if outdir is None:
        outdir = os.path.dirname(hap)

    # Check existence of input files
    paths_to_check = [hap, sample, ref_hap, ref_legend, recombination_map,
                      outdir]
    if to_impute is not None:
        paths_to_check += [to_impute]
    check_existence_of_paths(paths_to_check)

    if basename is None:
        basename = os.path.basename(hap).split(".gz")[0].split(".hap")[0]

    if suffix is None:
        suffix = ".impute2.{}_{}".format(from_bp, to_bp)

    imputed_hap = os.path.join(outdir, basename + suffix)
    cmd = [impute2_exe,
           "-use_prephased_g",
           "-known_haps_g", hap,
           "-sample_g",     sample,
           "-h",            ref_hap,
           "-l",            ref_legend,
           "-m",            recombination_map,
           "-Ne",           str(Ne),
           "-int",          str(from_bp), str(to_bp),
           "-buffer",       str(buffer_kb),
           "-o",            imputed_hap]

    if to_impute is not None:
        cmd += ["-include_snps", to_impute]

    if allow_large_regions:
        cmd += ["-allow_large_regions"]

    if chromosome in pysnpimpute.X_REGIONS:
        cmd += ["-chrX"]
        if chromosome.startswith("X_PAR"):
            cmd += ["-Xpar"]

    run_cmd_and_check(cmd, logger=logger)

    return imputed_hap
