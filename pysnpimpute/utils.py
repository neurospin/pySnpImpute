# -*- coding: utf-8 -*-

"""
A set of common utility functions.
"""

import os
import subprocess
import time
import argparse
import logging
import json
import glob

from collections import OrderedDict

import numpy
import pandas
import pysnpimpute

from pysnpimpute import (POSITIONS_OF_X_REGION_OF_BUILD,
                         CENTROMERE_POSITIONS_OF_CHROMOSOME_OF_BUILD)
from pysnpimpute.exceptions import (MissingSoftware,
                                    BadChromosomeName,
                                    BadGenomicBuildVersion,
                                    NotEnoughHaplotypes,
                                    AutosomeNotImputable)


def check_existence_of_paths(paths):
    """
    Raise an exception if any of the path in 'paths' (file or dir) does
    not exist.
    """
    for path in paths:
        if not os.path.exists(path):
            raise ValueError("File or directory does not exist: %s" % path)


def check_chromosome_name(chromosome):
    """
    Check that the chromosome name is handled.
    """
    if chromosome not in pysnpimpute.CHROMOSOMES:
        raise BadChromosomeName(chromosome)


def check_genomic_build_version(genomic_build):
    """
    Check that the passed <build_version> is handled.
    """
    if genomic_build not in pysnpimpute.GENOMIC_BUILDS:
        raise BadGenomicBuildVersion(genomic_build)


def extract_chromosome(bfile, chromosome, build=None, outdir=None,
                       basename=None, plink_exe="plink", logger=None):
    """
    Extract a chromosome or a X region from the dataset in Plink BED format.

    Parameters
    ----------
    chromosome: str
        Name of chromosome or X region to extract.
        Accepted names: "1", ..., "22", "X_PAR1", "X_nonPAR", "X_PAR2".
    bfile: str
        Path without extension to the dataset in PLINK BED format.
    build: str, default None
        Genomic build version of the <bfile> data (e.g. 'hg19'). Required if
        <chromosome> is a X region (the positions where to split depend of the
        genomic build).
    outdir: str, default None
        Path to directory where to output. By default in <bfile> directory.
    basename: str, default None
        Output filename without extension.
        By default input filename + ".<chromosome>".
    plink_exe: str, default "plink"
        Path to the Plink executable or alias if it's in $PATH.
    logger: logging object, defaut None
        To activate logging, pass a logging object.
    """

    # Check arguments and installation of required softwares
    check_chromosome_name(chromosome)
    if chromosome in pysnpimpute.X_REGIONS:
        check_genomic_build_version(build)
    check_installation_of_required_softwares(dict(Plink=plink_exe))

    # Output in <bfile> directory if <outdir> is not given.
    if outdir is None:
        outdir = os.path.dirname(bfile)

    # Check existence of input files
    paths_to_check = [bfile + ".bed", bfile + ".bim", bfile + ".fam", outdir]
    check_existence_of_paths(paths_to_check)

    if basename is None:
        basename = os.path.basename(bfile) + ".chr%s" % chromosome

    plink_chrom = "X" if chromosome in pysnpimpute.X_REGIONS else chromosome
    bfile_chrom = os.path.join(outdir, basename)
    cmd = [plink_exe,
           "--bfile",  bfile,
           "--chr",    plink_chrom,
           "--make-bed",
           "--out",    bfile_chrom,
           "--noweb"]

    # If the chromosome is a X region, add the positions where to cut
    if chromosome in pysnpimpute.X_REGIONS:
        from_bp, to_bp = POSITIONS_OF_X_REGION_OF_BUILD[build][chromosome]
        cmd += ["--from-bp", str(from_bp), "--to-bp",   str(to_bp)]

    # Run command with 'run_cmd' and not with 'run_cmd_and_check'
    # We check after the run if there are variants left
    run_cmd(cmd, logger=logger)

    # Determine the number of variants left: nb of lines of the .bim file
    path_bim = bfile_chrom + ".bim"
    if os.path.isfile(path_bim):
        with open(path_bim) as f:
            nb_variants_left = sum(1 for _ in f)
    else:
        nb_variants_left = 0
        bfile_chrom = None

    return bfile_chrom, nb_variants_left


def chunk_positions(region, known_hap_positions, ref_panel_positions,
                    chunksize_kb, N):
    """
    Take 2 lists of positions, known and ref panel variant positions, and
    compute intervals that define chunks.
    These chunks should:
        - not overlap
        - contain all haplotype positions
        - have at least <N> known haplotypes
        - have a minimum size of <chunksize_kb> in Kb (NOT IN BASEPAIR)

    region: name of the region being chunked. Used to generate clear errors.
            e.g. 'first arm of chr1', 'X_PAR1' etc
    """

    # Convert chunksize (Kb) to basepair
    chunksize = chunksize_kb * 1000

    known_hap_positions = numpy.asarray(known_hap_positions)
    ref_panel_positions = numpy.asarray(ref_panel_positions)

    nb_known_hap = len(known_hap_positions)
    if nb_known_hap < N:
        raise NotEnoughHaplotypes(region=region, nb_known=nb_known_hap,
                                  nb_required=N)

    # List of chunks, a chunk is a couple: (<start>, <end>)
    chunks = []
    ref_counts = []
    hap_counts = []

    while len(ref_panel_positions) > 0:

        start = ref_panel_positions[0]
        end = max(known_hap_positions[N-1], start + chunksize - 1)

        # Count the number of hap and ref hap in the chunk
        nb_hap = len(known_hap_positions[known_hap_positions <= end])
        nb_ref = len(ref_panel_positions[ref_panel_positions <= end])

        # Left to chunk
        ref_panel_positions = ref_panel_positions[ref_panel_positions > end]
        known_hap_positions = known_hap_positions[known_hap_positions > end]

        # If there is not enough variants left to create another chunk
        # merge with the current chunk
        nb_hap_left = len(known_hap_positions)
        nb_ref_left = len(ref_panel_positions)
        interval_left = ref_panel_positions[-1] - ref_panel_positions[0]
        if nb_hap_left < N or interval_left < chunksize:
            end = ref_panel_positions[-1]
            ref_panel_positions = []
            nb_hap += nb_hap_left
            nb_ref += nb_ref_left

        # Force type int, instead of numpy.int, otherwise by comparing them we
        # get numpy._bool and not bool, which is not handled by Hopla
        chunk = (int(start), int(end))
        chunks += [chunk]
        hap_counts += [nb_hap]
        ref_counts += [nb_ref]

    return chunks, hap_counts, ref_counts


def chromosome_segmentation(chromosome, known_hap, ref_legend, build,
                            chunksize_kb=5000, N=200, logger=None):
    """
    Compute intervals that define chromosome chunks to impute separately.
    These chunks should:
        - not overlap
        - contain all haplotype positions
        - not overlap with the centromere region of the chromosome
        - have at least <N> known haplotypes
        - have a minimum size of <chunksize_kb> Kb (NOT IN BASEPAIR)
    """

    # Check input arguments
    check_chromosome_name(chromosome)
    check_genomic_build_version(build)
    check_existence_of_paths([known_hap, ref_legend])

    # Load the positions of the known haplotypes
    df_hap = pandas.read_csv(known_hap, sep=" ", usecols=[2], header=None,
                             names=["pos"])

    # Load the positions of the reference haplotypes
    df_ref = pandas.read_csv(ref_legend, sep=" ", usecols=[1], header=0,
                             names=["pos"])

    # For the X regions X_PAR1 and X_PAR2 there are no centromere to handle
    # The chunking is simpler.
    if chromosome in {pysnpimpute.X_PAR1, pysnpimpute.X_PAR2}:
        chunks, known_hap_counts, ref_hap_counts = \
            chunk_positions(chromosome, df_hap["pos"].tolist(),
                            df_ref["pos"].tolist(), chunksize_kb, N)

    else:
        # For chromosome or regions where there is a centromere to handle,
        # chunking is done separately for both sides of the centromere.

        # Start and end positions of the centromere
        # chromosome structure: arm 1 | centromere | arm 2
        c_start, c_end = \
            CENTROMERE_POSITIONS_OF_CHROMOSOME_OF_BUILD[build][chromosome]

        # Dicts: Map <arm> -> list of haplotype positions
        known_hap_pos_of_arm = {
            "arm 1": df_hap[df_hap["pos"] < c_start]["pos"].tolist(),
            "arm 2": df_hap[df_hap["pos"] > c_end]["pos"].tolist()
        }

        ref_hap_pos_of_arm = {
            "arm 1": df_ref[df_ref["pos"] < c_start]["pos"].tolist(),
            "arm 2": df_ref[df_ref["pos"] > c_end]["pos"].tolist()
        }

        # Chunk both sides of the centromere separately

        # List to accumulate chunks and known/ref haplotype counts
        chunks, known_hap_counts, ref_hap_counts = [], [], []

        for arm in known_hap_pos_of_arm:
            region = "%s of chromosome %s" % (chromosome, arm)
            try:
                arm_chunks, arm_known_hap_counts, arm_ref_hap_counts = \
                    chunk_positions(region, known_hap_pos_of_arm[arm],
                                    ref_hap_pos_of_arm[arm], chunksize_kb, N)
                chunks.extend(arm_chunks)
                known_hap_counts.extend(arm_known_hap_counts)
                ref_hap_counts.extend(arm_ref_hap_counts)

            except NotEnoughHaplotypes as e:
                if logger is not None:
                    logger.error(e.message)

        # If none of the chromosome arms have enough haplotypes
        if len(chunks) == 0:
            raise AutosomeNotImputable(chromosome)

    return chunks, known_hap_counts, ref_hap_counts


def load_ref_panel_description(ref_panel):
    """
    Load the reference panel table as a Pandas DataFrame.

    Parameters
    ----------
    ref_panel: str
        Path to the table file describing the reference panel.
        It relates a chromosome name to 4 references files in Impute2 format:
         - reference haplotyples
         - reference legend
         - reference samples
         - recombination map
        The paths can be absolute or relative to the table file.
        File structure: one row per chromosome, tab-separated fields:
        <CHROM> <REF HAP> <REF LEGEND> <REF SAMPLE> <RECOMBINATION MAP>
        Only the following chromosome names are accepted:
        '1', ..., '22', 'X_PAR1', 'X_nonPAR' and 'X_PAR2'

    Return
    ------
    df_ref_panel: Pandas DataFrame, maps the 4 reference files for each chrom.
    """

    # Load the ref_panel txt file
    path_cols = ["ref_hap", "ref_legend", "ref_sample", "recombination_map"]
    df_ref_panel = pandas.read_csv(ref_panel, header=None, sep="\t", dtype=str,
                                   names=["chrom"] + path_cols, index_col=0)
    # Force index to be of type str
    df_ref_panel.index = df_ref_panel.index.astype(str)

    # If the paths are relative, replace with absolute paths
    ref_dir = os.path.dirname(os.path.abspath(ref_panel))
    to_abs_path = lambda p: p if os.path.isabs(p) else os.path.join(ref_dir, p)
    df_ref_panel[path_cols] = df_ref_panel[path_cols].applymap(to_abs_path)

    return df_ref_panel


def create_phased_dataset_description(phasing_outdir):
    """
    Create a "phased_dataset.txt" file in <phasing_outdir> mapping the paths
    to the phased data files to the corresponding chromosome names and sample
    files.
    It serves as input to the hopla_imputation.py script.
    """

    path_desc = os.path.join(phasing_outdir, "phased_dataset.txt")
    with open(path_desc, "w") as f:
        for c in pysnpimpute.ORDERED_CHROMOSOMES:
            regex = os.path.join(phasing_outdir, "chr%s/*.phased.hap.gz" % c)
            paths = glob.glob(regex)

            if len(paths) == 0:
                continue
            if len(paths) > 1:
                raise ValueError("Multiple files maching %s ." % regex)

            _, name_hap = os.path.split(paths[0])
            name_sample = name_hap[:-len(".hap.gz")] + ".sample"

            l = [c, "chr%s/%s" % (c, name_hap), "chr%s/%s" % (c, name_sample)]
            f.write("\t".join(l) + "\n")

    return path_desc


def load_phased_dataset_description(phased_dataset):
    """
    Load the phased dataset table file as a Pandas DataFrame.

    Parameters
    ----------
    phased_dataset: str
        Path the table file describing the phased data to impute.
        It relates a chromosome name to 2 files:
         - phased haplotypes to impute in Impute2 format
         - samples to impute in Impute2 format
        The paths can be absolute or relative to the table file.
        File structure: one row per chromosome, tab-separated fields:
        <CHROM> <HAP> <SAMPLE>
        Only the following chromosome names are accepted:
        '1', ..., '22', 'X_PAR1', 'X_nonPAR' and 'X_PAR2'
    """
    # Load the ref_panel txt file
    path_cols = ["hap", "sample"]
    df_dataset = pandas.read_csv(phased_dataset, header=None, sep="\t",
                                 dtype=str, names=["chrom"] + path_cols,
                                 index_col=0)
    # Force index to be of type str
    df_dataset.index = df_dataset.index.astype(str)

    # If the paths are relative, replace with absolute paths
    ref_dir = os.path.dirname(os.path.abspath(phased_dataset))
    to_abs_path = lambda p: p if os.path.isabs(p) else os.path.join(ref_dir, p)
    df_dataset[path_cols] = df_dataset[path_cols].applymap(to_abs_path)

    return df_dataset


def is_file(filepath):
    """ Check file's existence. Used as 'type' with the argparse module.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError("File does not exist: %s" % filepath)
    return filepath


def is_dir(dirpath):
    """ Check direcory's existence - argparse 'type' argument.
    """
    if not os.path.isdir(dirpath):
        raise argparse.ArgumentError("Directory does not exist: %s" % dirpath)
    return dirpath


def run_cmd(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, logger=None):
    """
    cmd: list of str
        Subprocess-like command, e.g. ["ls", "-l"]
    """
    # Command as a string for logging
    str_cmd = " ".join(cmd)

    if logger is not None:
        logger.info('cmd = "%s"' % str_cmd)

    start = time.strftime("%Y-%m-%d_%H:%M:%S")
    process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
    cout, cerr = process.communicate()
    exitcode = process.returncode
    end = time.strftime("%Y-%m-%d_%H:%M:%S")

    # Error message is stored as a list of lines
    stderr_lines = cerr.split("\n")

    couples = [("cmd", str_cmd), ("exitcode", exitcode),
               ("start", start), ("end", end), ("stderr", stderr_lines)]
    exec_metadata = OrderedDict(couples)

    if logger is not None:
        logger.info("exec_metadata = %s" % json.dumps(exec_metadata, indent=4))

    return exec_metadata


def run_cmd_and_check(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                      logger=None):
    """
    cmd: list of str
        Subprocess-like command, e.g. ["ls", "-l"]
    """
    exec_metadata = run_cmd(cmd, stdout=stdout, stderr=stderr, logger=logger)

    if exec_metadata["exitcode"] != 0:
        raise Exception("Command failed: {} .".format(" ".join(cmd)))


def check_installation_of_required_softwares(exe_of_software):
    """
    Parameters
    ----------
    exe_of_software: dict
        Map name of required software to the path or alias of the executable.
        An alias is enough if the executable is in the $PATH.

    Raise an exception if at least one the executable is not available.
    """
    for software_name in exe_of_software:
        exe_path_or_alias = exe_of_software[software_name]
        exec_metadata = run_cmd(["which", exe_path_or_alias])

        if exec_metadata["exitcode"] != 0:
            raise MissingSoftware(software_name, exe_path_or_alias)


def create_logger(prefix, log_dir=None, logger_name=None):
    """
    Create a logger with console logging (info level) + log file (debug level).
    If log_dir is None then logging directory is current directory.

    Log path is <log_dir>/<prefix>_%Y-%m-%d_%H:%M:%S.log
    """

    if logger_name is None:
        logger_name = __name__

    logger = logging.getLogger(logger_name)

    # Stop here if logger already has handlers
    if len(logger.handlers) != 0:
        return logger

    logger.setLevel(logging.DEBUG)

    # Log path
    log_filename = "%s_%s.log" % (prefix, time.strftime("%Y-%m-%d_%H:%M:%S"))
    if log_dir is not None:
        log_path = os.path.join(log_dir, log_filename)
    else:
        log_path = log_filename

    # File logger
    file_handler = logging.FileHandler(log_path)
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    # Console logger
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info("Path to log file: %s" % log_path)

    return logger
