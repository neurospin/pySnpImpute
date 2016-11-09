# -*- coding: utf-8 -*-

"""
Estimate the haplotypes prior to running the imputation using Shapeit2.
Phasing before imputing makes the imputation much faster.

The module requires that Plink and Shapeit2 are installed.
"""

import os

import pysnpimpute

from pysnpimpute.exceptions import NotEnoughGenotypes
from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_existence_of_paths,
                               extract_chromosome,
                               run_cmd_and_check)
from pysnpimpute.strand_alignment import align_dataset_to_ref_panel


def phase(chromosome, bfile, ref_hap, ref_legend, ref_sample,
          recombination_map, nb_cpus, N=200, build=None, Ne=20000, outdir=None,
          basename=None, suffix=".phased", plink_exe="plink",
          shapeit2_exe="shapeit", logger=None):
    """
    Align data of <chromosome> to reference panel (strand alignment) and
    phase. Both operations are done using Shapeit2 and the reference panel.

    Parameters
    ----------
    chromosome: str
        Name of chromosome or X region to phase.
        Accepted names: "1", ..., "22", "X_PAR1", "X_nonPAR", "X_PAR2".
    bfile: str
        Path without extension to the dataset in PLINK BED format.
    ref_hap, ref_legend, ref_sample: str
        Path to reference panel file in Shapeit2/Impute2 format.
    recombination_map: str
        Path to the recombination map required by Shapeit2 and Impute2.
    nb_cpus: int
        Number of CPUs to use in the phasing step (Shapeit2 --thread argument).
    N: int, default 200
        Minimum number of variants (i.e. known genotypes) to pursue with
        phasing after extracting the variants of the chromosome to phase.
        Otherwise raise NotEnoughGenotypes exception.
    build: str, default None
        Genomic build version of the <bfile> data (e.g. 'hg19'). Required if
        <chromosome> is a X region (the positions where to split depend of the
        genomic build).
    Ne: int, default 20000
        Shapeit2 '--effective-size' parameter. See also Impute2 'Ne' parameter.
    outdir: str, default None.
        Path to directory where to output. By default in <bfile> directory.
    plink_exe: str, default "plink"
        Path to the Plink executable or alias if it's in $PATH.
    shapeit2_exe: str, default "shapeit"
        Path to the shapeit2 executable or alias if it's in $PATH.
    logger: logging object, defaut None.
        To activate logging, pass a logging object.

    Raise
    -----
    NotEnoughGenotypes: raise if there is less than <N> variants left after
                        extracting the variants of the chromosome to phase.
    """

    # Check that Shapeit2 is installed
    check_installation_of_required_softwares(dict(Plink=plink_exe,
                                                  Shapeit2=shapeit2_exe))

    if outdir is None:
        outdir = os.path.dirname(bfile)

    if basename is None:
        basename = os.path.basename(bfile) + ".chr%s" % chromosome

    # Check existence of input files
    paths_to_check = [bfile + ".bed", bfile + ".bim", bfile + ".fam", outdir,
                      ref_hap, ref_legend, ref_sample, recombination_map]
    check_existence_of_paths(paths_to_check)

    # Extract chromosome (or X region) to phase for the <bfile> dataset.
    bfile_chrom, nb_variants = extract_chromosome(bfile=bfile,
                                                  chromosome=chromosome,
                                                  build=build,
                                                  outdir=outdir,
                                                  basename=basename,
                                                  plink_exe=plink_exe,
                                                  logger=logger)

    if nb_variants < N:
        raise NotEnoughGenotypes(region=chromosome, nb_known=nb_variants,
                                 nb_required=N)

    # Align dataset to reference panel (strand alignment)
    bfile_aligned = align_dataset_to_ref_panel(outdir=outdir,
                                               bfile=bfile_chrom,
                                               ref_hap=ref_hap,
                                               ref_legend=ref_legend,
                                               ref_sample=ref_sample,
                                               basename=basename,
                                               suffix=".refAligned",
                                               plink_exe=plink_exe,
                                               shapeit2_exe=shapeit2_exe,
                                               logger=logger)

    # Output paths
    basename = os.path.basename(bfile_aligned)
    out_basepath = os.path.join(outdir, basename + suffix)
    hap_phased = out_basepath + ".hap.gz"
    sample = out_basepath + ".sample"
    log_basepath = out_basepath + ".log"

    # Shapeit2 call
    cmd = [shapeit2_exe,
           "--input-bed",      bfile_aligned,
           "--input-ref",      ref_hap, ref_legend, ref_sample,
           "--input-map",      recombination_map,
           "--output-max",     hap_phased, sample,
           "--output-log",     log_basepath,
           "--effective-size", str(Ne),
           "--thread",         str(nb_cpus)]
    if chromosome == pysnpimpute.X_nonPAR:
        cmd += ["--chrX"]

    run_cmd_and_check(cmd, logger=logger)

    return hap_phased, sample
