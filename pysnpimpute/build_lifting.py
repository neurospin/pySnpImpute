# -*- coding: utf-8 -*-

"""
A set of functions to allow lifting from one genomic assembly build to
another (.e.g. from hg18 to hg19). The core computation is done by pyliftover
which is a Python implementation of the LiftOver software.
"""

import os

import pandas
import pyliftover

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_genomic_build_version,
                               run_cmd_and_check,
                               check_existence_of_paths)


def lift_build(bfile, from_build, to_build, outdir=None, basename=None,
               plink_exe="plink", logger=None):
    """
    Convert one genomic build into another (e.g. hg18 to hg19) using
    pyliftover.

    Returns
    -------
    bfile_lifted: str
        Path to output Plink BED data in the target build:
    path_unlifted: str
        Path to the list of variants that could not be converted to target
        build: <outdir>/<basename>.unlifted
    """

    # Check that Plink is installed
    check_installation_of_required_softwares(dict(Plink=plink_exe))

    if outdir is None:
        outdir = os.path.dirname(bfile)

    check_existence_of_paths([bfile + ".bed", bfile + ".bim", bfile + ".fam",
                              outdir])

    if basename is None:
        basename = os.path.basename(bfile)
    basename += ".%s" % to_build

    check_genomic_build_version(from_build)
    check_genomic_build_version(to_build)

    # Dict to adapt chromosome names for pyliftover
    plink_to_ucsc = {str(n): "chr%i" % n for n in range(1, 23)}
    plink_to_ucsc.update({'23': 'chrX', '24': 'chrY', '25': 'chrX',
                          '26': 'chrM'})

    # pyliftover object for lifting
    liftover = pyliftover.LiftOver(from_build, to_build)

    # Load Plink .bim file as a DataFrame
    df = pandas.read_csv(bfile + ".bim", header=None, sep="\t",
                         names=["chrom", "rs_id", "_1", "pos", "_2", "_3"],
                         dtype={"chrom": str})

    # Function to apply to every row of the Pandas DataFrame
    def lift_coordinate(row):
        plink_chrom = row["chrom"]
        ucsc_chrom = plink_to_ucsc[plink_chrom]
        try:
            return liftover.convert_coordinate(ucsc_chrom, row['pos'])[0][1]
        except:
            return -1

    # Apply lifting to all rows
    df["pos"] = df.apply(lift_coordinate, axis=1)

    # Identify variants for which the lifting failed
    unlifted_variants = df["rs_id"][df["pos"] < 0].tolist()

    # Create a file of variants to exclude for Plink --exclude
    path_unlifted = os.path.join(outdir, basename + ".unlifted_variants")
    with open(path_unlifted, "w") as f:
        f.write("\n".join(unlifted_variants) + "\n")

    # Remove unlifted variants from the DataFrame
    df = df[df["pos"] > 0]

    # Create a file for Plink --update-map to update positions
    path_pos_update = os.path.join(outdir, basename + ".new_positions")
    df.to_csv(path_pos_update, header=False, index=False, sep="\t",
              columns=["rs_id", "pos"])

    # Use Plink to remove unlifted variants and update positions
    bfile_lifted = os.path.join(outdir, basename)
    cmd = [plink_exe,
           "--bfile",      bfile,
           "--make-bed",
           "--exclude",    path_unlifted,
           "--update-map", path_pos_update,
           "--out",        bfile_lifted]
    run_cmd_and_check(cmd, logger=logger)

    return bfile_lifted, path_unlifted
