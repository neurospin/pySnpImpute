# -*- coding: utf-8 -*-

"""
A set of quality control functions. Used either as preliminary steps for
removing 'bad' samples/variants or to assess quality.
"""

import os

import pandas

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_existence_of_paths,
                               run_cmd_and_check)


def clean_data(bfile, use_imputed_sex=False, min_maf=0.01, min_vcallrate=0.95,
               min_scallrate=0.95, min_hwe=10e-6, outdir=None, basename=None,
               plink_exe="plink", logger=None):
    """
    Filter variants and subjects from the <bfile> dataset using Plink.

    Remove a variant if any of the following criteria is true:
        - minor allele frequency < <min_maf>
        - variant call rate < <min_vcallrate>
        - Hardy-Weinberg Equilibrium < <min_hwe>
        - if multiple variants have the same locus, they are excluded

    Remove a subject if any of the following criteria is true:
        - subject/sample genotype call rate < <min_scallrate>
        - sex is missing (i.e. 0)
        - imputed sex differ from reported sex, if using reported sex

    Outpath: <outdir>/<basename> + .bed/bim/fam

    Return
    ------
    bfile_cleaned: str
        Path without extension to the output files (.bed/.bim/.fam).
        None if there are not variants left after filtering.
    nb_variants_left: int
        Number of variants left after filtering.
    """

    # Check that Plink is installed
    check_installation_of_required_softwares(dict(Plink=plink_exe))

    if bfile.endswith(".bed"):
        bfile = bfile[:-len(".bed")]

    # Check existence of input files
    paths_to_check = [bfile + ".bed", bfile + ".bim", bfile + ".fam"]
    check_existence_of_paths(paths_to_check)

    if outdir is None:
        outdir = os.path.dirname(bfile)

    if basename is None:
        basename = os.path.basename(bfile)

    out_basepath = os.path.join(outdir, basename)

    if use_imputed_sex:
        out_basepath += ".imputedsex"
        cmd = [plink_exe,
               "--bfile", bfile,
               "--impute-sex",
               "--make-bed",
               "--out",   out_basepath]
    else:
        cmd = [plink_exe,
               "--bfile", bfile,
               "--check-sex",
               "--out",   out_basepath]

    run_cmd_and_check(cmd, logger=logger)

    # Read Plink sex report and create a list of samples to remove
    path_report = out_basepath + ".sexcheck"
    df_sex = pandas.read_csv(path_report, delim_whitespace=True, dtype=str)

    if use_imputed_sex:
        df_sex = df_sex[df_sex["SNPSEX"] == "0"]
    else:
        df_sex = df_sex[df_sex["STATUS"] == "PROBLEM"]

    path_subjects_to_exclude = out_basepath + ".excluded_by_sex_qc"
    df_sex.to_csv(path_subjects_to_exclude, index=False, header=False,
                  sep="\t", columns=["FID", "IID"])

    bfile_to_qc = out_basepath if use_imputed_sex else bfile

    # Identify duplicates, i.e. variants that share the same locus
    cmd = [plink_exe,
           "--bfile", bfile_to_qc,
           "--list-duplicate-vars",
           "--out", out_basepath]
    run_cmd_and_check(cmd, logger=logger)

    # Filter the dataset
    bfile_cleaned = out_basepath + ".qc"
    path_duplicated_variants = out_basepath + ".dupvar"
    cmd = [plink_exe,
           "--bfile",   bfile_to_qc,
           "--remove",  path_subjects_to_exclude,
           "--exclude", path_duplicated_variants,
           "--maf",     str(min_maf),
           "--geno",    str(1-min_vcallrate),
           "--mind",    str(1-min_scallrate),
           "--hwe",     str(min_hwe),
           "--make-bed",
           "--out",     bfile_cleaned,
           "--noweb"]
    run_cmd_and_check(cmd, logger=logger)

    return bfile_cleaned


def mds_analysis(bfile, mind=0.05, mds_plot=4, outdir=None, plink_exe="plink",
                 logger=None):
    """
    Run MDS analysis using Plink. Set mds_plot for the number of directions.
    """

    # Check that Plink is installed
    check_installation_of_required_softwares(dict(Plink=plink_exe))

    if bfile.endswith(".bed"):
        bfile = bfile[:-len(".bed")]

    # Check existence of input files
    required_files = [bfile + ".bed", bfile + ".bim", bfile + ".fam"]
    check_existence_of_paths(required_files)

    if outdir is None:
        mdsfile = bfile + ".mds"
    else:
        mdsfile = os.path.join(outdir, os.path.basename(bfile) + ".mds")

    cmd = [plink_exe,
           "--bfile",    bfile,
           "--cluster",
           "--mind",     str(mind),
           "--mds-plot", str(mds_plot),
           "--out",      mdsfile,
           "--noweb"]
    run_cmd_and_check(cmd, logger=logger)

    return mdsfile
