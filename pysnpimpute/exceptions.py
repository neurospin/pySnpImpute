# -*- coding: utf-8 -*-

"""
Package exceptions, all deriving the base Exception: PySnpImputeError
"""

import pysnpimpute


class PySnpImputeError(Exception):
    """ Base exception type for the package.
    """
    def __init__(self, message):
        super(PySnpImputeError, self).__init__(message)


class MissingSoftware(PySnpImputeError):
    """
    Exception thrown when a required program is not installed or detected.
    """
    def __init__(self, software_name, exe_path_or_alias):
        message = ("Required software '%s' is not available as the setted "
                   "alias/path: %s" % (software_name, exe_path_or_alias))
        super(MissingSoftware, self).__init__(message)


class BadChromosomeName(PySnpImputeError):
    """
    Exception thrown when an unaccepted chromosome name is detected.
    """
    def __init__(self, chromosome):
        message = ("Bad chromosome name: {}, accepted values: {}"
                   .format(chromosome, pysnpimpute.CHROMOSOMES))
        super(BadChromosomeName, self).__init__(message)


class BadGenomicBuildVersion(PySnpImputeError):
    """
    Exception thrown when an unaccepted genomic build version is detected.
    """
    def __init__(self, genomic_build):
        message = ("Bad genomic build version: {}, supported: {}"
                   .format(genomic_build, pysnpimpute.GENOMIC_BUILDS))
        super(BadGenomicBuildVersion, self).__init__(message)


class NotEnoughGenotypes(PySnpImputeError):
    """
    Exception thrown when there are not enough measured genotypes in a region
    for phasing. A region can be a chromosome arm, a X region (e.g. X_PAR1),
    a chromosome chunk etc.
    """
    def __init__(self, region, nb_known, nb_required):
        message = ("Region {} cannot be phased. There are not enough "
                   "measured genotypes. Known: {}, required: {}"
                   .format(region, nb_known, nb_required))
        super(NotEnoughGenotypes, self).__init__(message)


class NotEnoughHaplotypes(PySnpImputeError):
    """
    Exception thrown when there are not enough known haplotypes in a region for
    imputation. A region can be a chromosome arm, a X region (e.g. X_PAR1),
    a chromosome chunk etc.
    """
    def __init__(self, region, nb_known, nb_required):
        message = ("Region {} was not imputed. There are not enough known "
                   "haplotypes. Known: {}, required: {}"
                   .format(region, nb_known, nb_required))
        super(NotEnoughHaplotypes, self).__init__(message)


class AutosomeNotImputable(PySnpImputeError):
    """
    Exception thrown when an autosome chromosome is not imputable at all
    because none of the chromosome arms have enough known haplotypes to be
    imputed.
    """
    def __init__(self, chromosome):
        message = ("Chromosome {} was not imputed at all. None of the "
                   "chromosome arms had enough known haplotypes to be imputed."
                   .format(chromosome))
        super(AutosomeNotImputable, self).__init__(message)
