# -*- coding: utf-8 -*-

"""
This module defines a pipeline to preprocess the dataset before it can be used
for phasing.
"""

import os

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_existence_of_paths)
from pysnpimpute.qc import clean_data
from pysnpimpute.build_lifting import lift_build


def preprocess(bfile, outdir=None, from_build=None, to_build=None,
               use_imputed_sex=False, min_maf=0.01, min_vcallrate=0.95,
               min_scallrate=0.95, min_hwe=10e-6, basename=None,
               plink_exe="plink", logger=None):
    """
    Preprocess the <bfile> dataset.

    Steps:
    - Lift the data to the desired genomic build, if requested. Required if
      the dataset to phase/impute is not already in the same build as the
      reference panel. Done using pyliftover.
    - Replace genders with imputed genders if requested (Plink --impute-sex)
    - Remove duplicated variants (Plink --list-duplicate-vars and --exclude)
    - Remove variants for which any of the following criteria is true:
        - minor allele frequency < <min_maf>
        - variant call rate < <min_vcallrate>
        - Hardy-Weinberg Equilibrium < <min_hwe>
    - Remove samples for which any of the following criteria is true:
        - sample genotype call rate < <min_scallrate>
        - sex is missing (i.e. 0)
        - imputed sex differ from reported sex, if using reported sex

    Parameters
    ----------
    bfile: str
        Path without extension to the dataset in PLINK BED format.
    outdir: str, default None.
        Path to directory where to output. By default in <bfile> directory.
    from_build, to_build: str, default None
        To change the build version (e.g. from hg18 to hg19), set the source
        and target builds. For imputation, if the data to impute is not
        already in the same build as the reference panel, a lift is required.
    use_imputed_sex: bool, default False
        If True replace reported sex with imputed sex (imputed with Plink).
    min_maf: float, default 0.01
        Minimum minor allele frequency for a variant to be kept.
    min_vcallrate: float, default 0.95
        Minimum genotyping rate for a variant to be kept.
    min_scallrate: float, default 0.95
        Minimum genotyping rate for a subject/sample to be kept.
    min_hwe: float, default 10e-6
        Minimum Hardy-Weinberg Equilibrium value for a SNP to be kept.
        Path to the liftOver executable or alias if it's in $PATH.
    plink_exe: str, default "plink"
        Path to the Plink executable or alias if it's in $PATH.
    logger: logging object, defaut None.
        To activate logging, pass a logging object.
    """

    check_installation_of_required_softwares(dict(Plink=plink_exe))

    # Output in <bfile> directory if <outdir> is not given.
    if outdir is None:
        outdir = os.path.dirname(bfile)

    if basename is None:
        basename = os.path.basename(bfile)

    # Create outdir if it does not exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Check existence of paths
    paths_to_check = [outdir, bfile + ".bed", bfile + ".bim", bfile + ".fam"]
    check_existence_of_paths(paths_to_check)

    # STEP 1 - Convert build if requested
    if (from_build, to_build) is not (None, None):
        bfile_to_qc, _ = lift_build(bfile=bfile,
                                    from_build=from_build,
                                    to_build=to_build,
                                    outdir=outdir,
                                    basename=basename,
                                    plink_exe=plink_exe,
                                    logger=logger)
    else:
        bfile_to_qc = bfile

    # STEP 2 - Quality control
    basename = os.path.basename(bfile_to_qc)
    bfile_qc = clean_data(bfile=bfile_to_qc,
                          use_imputed_sex=use_imputed_sex,
                          min_maf=min_maf,
                          min_vcallrate=min_vcallrate,
                          min_scallrate=min_scallrate,
                          min_hwe=min_hwe,
                          outdir=outdir,
                          basename=basename,
                          plink_exe=plink_exe,
                          logger=logger)

    return bfile_qc
