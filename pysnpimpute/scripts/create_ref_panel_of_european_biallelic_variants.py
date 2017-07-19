# -*- coding: utf-8 -*-

DOC = """
Download the 1000 Genomes phase 3 reference panel in VCF format, along with
recombination maps in Impute2 format, filter all the VCF files and convert to
Impute2 format, keeping only variants that are bi-allelic, of type SNP or
INDEL and that have a minor allele count >= 2 in the EUR super-population.

The idea is to have a smaller reference panel to decrease the cost of
imputation.

Requirements:
- The filtering is done using Bcftools, it should be installed.

Note:
- The X chromosome is splitted in 3 regions: X_PAR1, X_nonPAR, X_PAR2
"""

import os
import argparse
import multiprocessing

import pandas
import pysnpimpute

from pysnpimpute.utils import (check_installation_of_required_softwares,
                               check_existence_of_paths,
                               create_logger,
                               run_cmd_and_check)
from pysnpimpute import POSITIONS_OF_X_REGION_OF_BUILD


def download_1000g_phase3_ref_panel(download_dir):
    """
    Download the 1000 Genomes Project phase 3 reference for the autosomes and
    the X chromosome in <download_dir>.
    """

    # Create <download_dir> if it does not exist
    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    # wget download options
    #   - do not redownload if already downloaded
    #   - if it fails/stops, retry and continue previous download
    wget_options = ["-c", "-N", "-t", "10"]

    base_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

    # List all the URLs to download
    urls = []

    # Sample panel url
    filename = "integrated_call_samples_v3.20130502.ALL.panel"
    urls += [os.path.join(base_url, filename)]
    path_sample = os.path.join(download_dir, filename)

    # Dict mapping <chromosome> -> <path vcf>
    vcf_of_chrom = dict()

    # Add autosomes url
    for chrom in pysnpimpute.ORDERED_AUTOSOMES:
        filename = ("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a"
                    ".20130502.genotypes.vcf.gz" % chrom)
        url = os.path.join(base_url, filename)
        urls += [url, url + ".tbi"]
        vcf_of_chrom[chrom] = os.path.join(download_dir, filename)

    # X chromosome url
    filename = ("ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b"
                ".20130502.genotypes.vcf.gz")
    url = os.path.join(base_url, filename)
    urls += [url, url + ".tbi"]
    vcf_of_chrom["X"] = os.path.join(download_dir, filename)

    # Write vcf_urls.txt with all urls to download
    path_urls = os.path.join(download_dir, "vcf_urls.txt")
    with open(path_urls, "w") as f:
        f.write("\n".join(urls) + "\n")

    # Run downloading
    cmd = ["wget", "-i", path_urls, "-P", download_dir] + wget_options
    run_cmd_and_check(cmd)

    return path_sample, vcf_of_chrom


def download_recombination_maps(download_dir):
    """
    Download recombination maps for build hg19 (b37) in Impute2 format.
    """

    # Create <download_dir> if it does not exist
    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    # wget download options
    #   - do not redownload if already downloaded
    #   - if it fails/stops, retry and continue previous download
    wget_options = ["-c", "-N", "-t", "3"]

    base_url = "http://mathgen.stats.ox.ac.uk/impute/"
    filename = "genetic_maps_b37.tgz"
    url = os.path.join(base_url, filename)
    path = os.path.join(download_dir, filename)
    cmd = ["wget", url, "-P", download_dir] + wget_options
    run_cmd_and_check(cmd)

    # Extract files from compressed archive
    cmd = ["tar", "xzf", path, "--strip-components=1", "-C", download_dir]
    run_cmd_and_check(cmd)

    # Dict mapping <chrom> -> <path recombination map>
    rmap_of_chrom = dict()
    for chrom in pysnpimpute.CHROMOSOMES:
        fname = "genetic_map_chr%s_combined_b37.txt" % chrom
        rmap_of_chrom[chrom] = os.path.join(download_dir, fname)
    check_existence_of_paths(rmap_of_chrom.values())

    return rmap_of_chrom


def list_european_samples(path_panel, path_european_samplist):
    """
    Create a file listing the IDs of "EUR" samples.
    """
    df = pandas.read_csv(path_panel, sep="\t", usecols=["sample", "super_pop"])
    df = df[df["super_pop"] == "EUR"]  # Remove non european
    df.to_csv(path_european_samplist, index=False, header=False,
              columns=["sample"])

    return path_european_samplist


def filter_vcf_and_convert_to_impute2(chrom, path_vcf, path_european_samplist,
                                      outdir, bcftools_exe="bcftools"):
    """

    Parameters
    ----------
    chrom: str
        Name of chrom or X region to filter.
        Accepted names: "1", ..., "22", "X_PAR1", "X_nonPAR", "X_PAR2".
    path_vcf: str
        Path to the .vcf.gz of the reference panel to filter and convert.
    path_european_samplist: str
        Path to the list of EUR subjects.
    outdir: str, default None
        Path to directory where to output.
    bcftools: str, default "bcftools"
        Path to the 'bcftools' executable or alias if it's in $PATH.
    """

    # Check that bcftools is installed.
    check_installation_of_required_softwares(dict(bcftools=bcftools_exe))

    basename = "1000G_P3.ALLsamples.biallelic_EUR_variants.chr%s" % chrom
    out_basepath = os.path.join(outdir, basename)

    # Filter the .vcf.gz to keep only EUR subjects, to identify the variants
    # that are bi-allelic and found at least twice in this population (Minor
    # allele count >= 2).
    path_vcf_filtered = out_basepath + ".temp.vcf"
    cmd = [bcftools_exe,
           "view",           path_vcf,
           "--samples-file", path_european_samplist,
           "--min-ac",       "2",
           "--max-alleles",  "2",
           "--types",        "snps,indels",
           "--output-file",  path_vcf_filtered,
           "--output-type",  "v"]

    if chrom in pysnpimpute.X_REGIONS:
        from_bp, to_bp = POSITIONS_OF_X_REGION_OF_BUILD["hg19"][chrom]
        cmd += ["--regions", "X:%i-%i" % (from_bp, to_bp)]
    run_cmd_and_check(cmd)

    # Convert the ref .vcf.gz to Impute2 and keep only the variants from the
    # filtered .vcf.gz. Keep rs ids if possible, instead of CHROM_POS_REF_ALT
    cmd = [bcftools_exe,
           "convert",           path_vcf,
           "--regions-file",    path_vcf_filtered,
           "--haplegendsample", out_basepath]
    run_cmd_and_check(cmd)

    extensions = [".hap.gz", ".legend.gz", ".samples"]
    rhap, rlegend, rsample = map(lambda x: out_basepath + x, extensions)

    # Remove intermediate filtered VCF file
    os.remove(path_vcf_filtered)

    # Remove useless '.samples' file (not well formed by default)
    os.remove(rsample)

    return chrom, rhap, rlegend


def unpacker(kwargs):
    return filter_vcf_and_convert_to_impute2(**kwargs)


def main(outdir, nb_processes, bcftools_exe):
    """
    Download the 1000 Genomes phase 3 reference panel, along with recombination
    maps in Impute2 format, filter all the VCF files and convert to Impute2
    format, keeping only variants that are bi-allelic, of type SNP or INDEL
    and that have a minor allele count >= 2 in the EUR super-population.
    """

    # Create outdir if it does not exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    logger = create_logger(prefix="create_ref_panel", log_dir=outdir)

    # Download the 1000 Genomes Project reference panel files
    vcf_dir = os.path.join(outdir, "VCF")
    logger.info("Downloading 1000 Genomes Project phase 3 reference panel "
                "to %s . ~17GO to download, might take a while..." % vcf_dir)
    path_panel, vcf_of_chrom = download_1000g_phase3_ref_panel(vcf_dir)

    # Download the recombination maps
    impute2_dir = os.path.join(outdir, "Impute2_european_variants")
    logger.info("Downloading recombination maps in Impute2 format to %s . "
                "~50MO to donwload." % impute2_dir)
    rmap_of_chrom = download_recombination_maps(impute2_dir)

    # Convert the reference panel sample file in Shapeit2/Impute2 format
    # It implies replacing "male"/"female" by "1"/"2"
    df_panel = pandas.read_csv(path_panel, sep="\t", usecols=[0, 1, 2, 3])
    sex_map = {"male": "1", "female": "2"}
    df_panel.gender = df_panel.gender.apply(lambda x: sex_map[x])
    path_panel_impute2 = os.path.join(impute2_dir, "1000G_P3.ALLsamples.panel")
    df_panel.to_csv(path_panel_impute2, index=False, sep="\t")

    # Create a list of EUR subjects
    path_european_samplist = os.path.join(impute2_dir, "european.samples")
    list_european_samples(path_panel_impute2, path_european_samplist)

    # List to accumulate kwargs for parallel processing
    list_kwargs = []
    for chrom in pysnpimpute.CHROMOSOMES:
        # The X chromosome is in one VCF file but will be converted to
        # Impute2 format in 3 files: X_PAR1, X_nonPAR, X_PAR2
        vcf_chrom = "X" if chrom in pysnpimpute.X_REGIONS else chrom
        kwargs = dict(chrom=chrom, path_vcf=vcf_of_chrom[vcf_chrom],
                      path_european_samplist=path_european_samplist,
                      outdir=impute2_dir)
        list_kwargs += [kwargs]

    # Run concurrently
    logger.info("Filtering VCF files and convert to Impute2 format. Keeping "
                "only variants that are bi-allelic, of type SNP or INDEL and "
                "that have a minor allele count >= 2 in the EUR "
                "super-population. Might take a while...")
    pool = multiprocessing.Pool(nb_processes)
    results = pool.map(unpacker, list_kwargs)
    pool.join()
    pool.close()

    # Create a table file describing the reference panel
    ref_panel_txt = os.path.join(impute2_dir, "ref_panel.txt")
    with open(ref_panel_txt, "w") as f:
        for chrom, rhap, rlegend in results:
            paths = [rhap, rlegend, path_panel_impute2, rmap_of_chrom[chrom]]
            filenames = map(lambda x: os.path.basename(x), paths)
            fields = [chrom] + filenames
            f.write("\t".join(fields) + "\n")

    return ref_panel_txt


def get_cmd_line_args():
    """
    Parse the command line arguments and return a dict mapping
    <argument name> -> <value>.
    """

    parser = argparse.ArgumentParser(description=DOC)

    dhelp = ("Path to directory where to download and convert the 1000 "
             "Genomes phase3 reference panel.")
    parser.add_argument("-o", "--outdir", required=True, metavar="<path>",
                        help=dhelp)

    nhelp = ("Number of processes to run in parallel. Each process filters "
             "and converts one chromosome of the reference panel.")
    parser.add_argument("-n", "--nb-processes", type=int, required=True,
                        metavar="<int>", help=nhelp)

    # Optional arguments
    ehelp = "Path to the 'bcftools' executable or alias if it's in $PATH."
    parser.add_argument("-e", "--bcftools-exe", default="bcftools",
                        metavar="<exe>", help=ehelp)

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)

    return kwargs


if __name__ == "__main__":

    kwargs = get_cmd_line_args()
    main(**kwargs)
