
Preprocessing
=============

The dataset to impute has to be preprocessed before it can be phased.

The script takes less than a minute to run.

Script steps:
    - Lift the data to the desired genomic build, if requested.
      Required if the dataset to phase/impute is not already in the same
      build as the reference panel. Done using pyliftover.
    - Replace genders with imputed genders if requested (using Plink ``--impute-sex``)
    - Remove duplicated variants (using Plink ``--list-duplicate-vars`` and ``--exclude``)
    - Remove variants for which any of the following criteria is true:
        - minor allele frequency < <min_maf>
        - variant call rate < <min_vcallrate>
        - Hardy-Weinberg Equilibrium < <min_hwe>
    - Remove samples (subjects) for which any of the following criteria is true:
        - sample call rate < <min_scallrate>
        - sex is unknown (i.e. 0)
        - SNP sex differ from reported sex, if not using imputed sex

Note: the lifting, if requested, is done before removing duplicated variants.
In fact 2 variants could have different loci in the raw data build and then
be merged in a more recent build (target build) and thus become duplicates.
To make sure these cases don't appear, we have to run detection of duplicates
after lifting the build.

=============================== ======================================================================== ==========
 *preprocess.py*
-------------------------------------------------------------------------------------------------------------------
          Parameter                                           Description                                 Default
=============================== ======================================================================== ==========
``-i, --bfile``                 Path without extension to the dataset in Plink bed/bim/fam format.       *required*
``-o, --outdir <path>``         Path to directory where to output.                                       *required*
``-l, --lift <from> <to>``      To change the build version (e.g. from hg18 to hg19), set the source     *None*
                                and target builds. If the data to phase/impute is not already in the
                                same build as the reference panel, lifting is required.
``-p, --use-imputed-sex``       Replace reported sex with imputed sex (imputed with*Plink --impute-sex). *False*
``-a, --min-maf <float>``       Minimum minor allele frequency.                                          *0.05*
``-k, --min-vcallrate <float>`` Minimum genotyping rate for a variant to be kept.                        *0.95*
``-m, --min-scallrate <float>`` Minimum genotyping rate for a sample to be kept.                         *0.95*
``-w, --min-hwe <float>``       Minimum Hardy-Weinberg equilibrium.                                      *10e-6*
``-e, --plink-exe <exe>``       Path to the Plink executable or alias if it's in $PATH.                  *'plink'*
``-z, --no-verbose``            Do not log.                                                              *False*
=============================== ======================================================================== ==========
