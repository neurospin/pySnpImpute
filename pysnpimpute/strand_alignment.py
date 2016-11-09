# -*- coding: utf-8 -*-

"""
A set of functions to check and correct strand alignment. If a variant is not
alignable it is removed. Alignment is required before running imputation.

The module requires Shapeit2 and Plink to be installed:
- Shapeit2 is used to check the strand alignment
- Plink is used to flip SNPs that are not aligned and to remove SNPs that are
  not alignable.
"""

import os
import tempfile
import shutil

import pandas

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_existence_of_paths,
                               run_cmd,
                               run_cmd_and_check)


def check_alignment_to_ref_panel(bfile, ref_hap, ref_legend, ref_sample,
                                 outdir=None, basename=None,
                                 suffix=".misaligned", shapeit2_exe="shapeit",
                                 logger=None):
    """
    Create a file listing the SNPs which strand is not aligned to reference
    panel using the "-check" option of Shapeit2.

    Parameters
    ----------
    bfile: str
        Path to input Plink binary file.
    ref_hap, ref_legend, ref_sample: str
        Path to reference panel file.
    outdir: str, default None
        Path to directory where to output. By default <bfile> directory.
    basename: str, default None
        Output basename, by default <bfile> basename.
    suffix: str, default ".misaligned"
        Suffix concatenated to <basename>.
    shapeit2_exe: str, default "shapeit"
        Path to the shapeit2 executable or alias if it's in $PATH.
    logger: logging object, defaut None.
        To activate logging, pass a logging object.

    Return
    ------
    snpfile_misaligned: str
        Path to snp file, <outdir>/<basename><suffix>
    """

    # Check that Shapeit2 is installed
    check_installation_of_required_softwares(dict(Shapeit2=shapeit2_exe))

    if outdir is None:
        outdir = os.path.dirname(bfile)

    if basename is None:
        basename = os.path.basename(bfile)

    # Create a temporary directory where shapeit2 can write the logs
    tmp_dir = tempfile.mkdtemp(prefix="checkAlignmentToRefPanel_")

    out_basepath = os.path.join(tmp_dir, basename)
    cmd = [shapeit2_exe,
           "-check",
           "--input-bed",  bfile,
           "--input-ref",  ref_hap, ref_legend, ref_sample,
           "--output-log", out_basepath]
    run_cmd(cmd, logger=logger)

    shapeit2_misaligned_log = out_basepath + ".snp.strand"
    if os.path.isfile(shapeit2_misaligned_log):
        df = pandas.read_csv(shapeit2_misaligned_log, sep="\t")
        misaligned_snps = list(set(df["main_id"]))
    else:
        misaligned_snps = []

    # Remove temporary directory
    shutil.rmtree(tmp_dir)

    snpfile_misaligned = os.path.join(outdir, basename + suffix)
    with open(snpfile_misaligned, "w") as f:
        f.write("\n".join(misaligned_snps) + "\n")

    return snpfile_misaligned


def align_snps(bfile, flip_snps, exclude_snps=None, outdir=None,
               basename=None, suffix=".aligned", plink_exe="plink",
               logger=None):
    """
    Align SNPs by flipping SNPs listed in the <flip_snps> file and excluding
    SNPs listed in <exclude_snps> using Plink.
    """

    # Check that Plink is installed
    check_installation_of_required_softwares(dict(Plink=plink_exe))

    if outdir is None:
        outdir = os.path.dirname(bfile)

    if basename is None:
        basename = os.path.basename(bfile)

    bfile_aligned = os.path.join(outdir, basename + suffix)
    cmd = [plink_exe,
           "--bfile",    bfile,
           "--flip",     flip_snps,
           "--make-bed",
           "--out",      bfile_aligned,
           "--noweb"]

    if exclude_snps is not None:
        cmd += ["--exclude", exclude_snps]

    run_cmd_and_check(cmd, logger=logger)

    return bfile_aligned


def align_dataset_to_ref_panel(outdir, bfile, ref_hap, ref_legend, ref_sample,
                               basename=None, suffix=".refAligned",
                               plink_exe="plink", shapeit2_exe="shapeit",
                               logger=None):
    """
    Align GWAS data (Plink format) of a chromosome to reference panel prior
    to phasing.
    The GWAS dataset should be in the same build as the reference panel.

    Steps:
        - detect alignment errors (i.e. SNPs not aligned) using shapeit2
        - flip all SNPs that are not aligned using Plink
        - detect persistent alignment errors (i.e. SNPs still not aligned)
          using shapeit2
        - remove these SNPs from the dataset using Plink

    Parameters
    ----------
    outdir: str
        Path to directory where to output.
    bfile: str
        Path to the GWAS data (Plink binary format) to align.
    ref_hap, ref_legend, ref_sample: str
        Path to reference panel file in shape2/impute2 format.
    logger: logging object, defaut None.
        To activate logging, pass a logging object.
    suffix: str
        Output filenames: <outdir>/<bfile name><suffix>.bed/bim/fam
    plink_exe: str, default "plink"
        Path to the Plink executable or alias if it's in $PATH.
    shapeit2_exe: str, default "shapeit"
        Path to the shapeit2 executable or alias if it's in $PATH.
    """

    # Check that Plink and Shapeit2 are installed
    check_installation_of_required_softwares(dict(Plink=plink_exe,
                                                  Shapeit2=shapeit2_exe))

    # Remove .bed extension if used
    if bfile.endswith(".bed"):
        bfile = bfile[:-len(".bed")]

    # Check input files and outdir existence
    required_paths = [outdir, bfile + ".bed", bfile + ".bim", bfile + ".fam",
                      ref_hap, ref_legend, ref_sample]
    check_existence_of_paths(required_paths)

    # Step 1 - Detect alignment errors
    misaligned_snps = check_alignment_to_ref_panel(bfile=bfile,
                                                   ref_hap=ref_hap,
                                                   ref_legend=ref_legend,
                                                   ref_sample=ref_sample,
                                                   basename=basename,
                                                   suffix=".misaligned",
                                                   shapeit2_exe=shapeit2_exe,
                                                   logger=logger)

    # STEP 2 - Flip all SNPs that are not aligned
    bfile_flipped = align_snps(bfile,
                               flip_snps=misaligned_snps,
                               outdir=outdir,
                               basename=basename,
                               suffix=suffix,
                               plink_exe=plink_exe,
                               logger=logger)

    # STEP 3 - Repeat alignment check on flipped dataset
    exclude_snps = check_alignment_to_ref_panel(bfile=bfile_flipped,
                                                ref_hap=ref_hap,
                                                ref_legend=ref_legend,
                                                ref_sample=ref_sample,
                                                basename=basename,
                                                suffix=".excluded_unalignable",
                                                shapeit2_exe=shapeit2_exe,
                                                logger=logger)

    # STEP 4 - Exclude SNPs still not aligned
    bfile_aligned = align_snps(bfile,
                               flip_snps=misaligned_snps,
                               exclude_snps=exclude_snps,
                               outdir=outdir,
                               basename=basename,
                               suffix=suffix,
                               plink_exe=plink_exe,
                               logger=logger)

    return bfile_aligned
