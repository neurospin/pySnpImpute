.. pySnpImpute documentation master file, created by
   sphinx-quickstart on Mon Jul  4 17:30:39 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pySnpImpute's documentation!
=======================================

*pySnpImpute* is a Python package meant to help with running genotype
imputation using a reference panel. For that matter, it defines
ready-to-run command line scripts for common task:
creating a reference panel, preprocessing, phasing and imputing.
*pySnpImpute* is essentially a wrapper to a set of standard tools:
*Bcftools*, *Plink*, *pyliftover*, *Shapeit2* and *Impute2*.


Installation
############
System or user installation::

  pip install pysnpimpute
  pip install --user pysnpimpute

Python dependencies that *pip* will automatically update/install if necessary:
  * *Numpy*
  * *Pandas*
  * *pyliftover*: a Python implementation of LiftOver
  * *Hopla*: Python package used to run multiple instance of a script,
    used to scale phasing and imputation

The binaries that the package requires are not distributed automatically,
you should install them, if not already the case. Each script will check
whether the tools it requires are installed and will throw an error if it
cannot find them. The scripts do not require the binaries to be in the
$PATH as you can optionnally pass the paths to the binaries as arguments.


==========  ========================================================== ==============================================================
   Tool                                Link                                                     Used for
==========  ========================================================== ==============================================================
*Bcftools*  `<https://samtools.github.io/bcftools>`_                   Filtering VCF files and conversion to Shapeit2/Impute2 format.
*Plink2*    `<http://www.cog-genomics.org/plink2>`_                    QC/filtering of input data and sex imputation.
*Shapeit2*  `<http://www.shapeit.fr>`_                                 Checking of strand alignment and estimation of the haplotypes.
*Impute2*   `<https://mathgen.stats.ox.ac.uk/impute/impute_v2.html>`_  Imputation of haplotypes.
==========  ========================================================== ==============================================================


Imputation steps
################

.. toctree::

  reference_panel
  dataset_to_impute
  preprocessing
  phasing
  imputation
  example_ADNI



.. Contents:

.. toctree::
   :maxdepth: 2



