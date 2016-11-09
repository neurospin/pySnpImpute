
Phasing
=======

The phasing, or estimation of haplotypes, is done with the help of Shapeit2.

The package provides 2 scripts for that matter:
  - *phase.py*: to run phasing for one chromosome only (or X region).
                It extracts the data of the selected chromosome (or X region)
                from the Plink preprocessed dataset, aligns the strands to
                the reference panel and phases the chosen chromosome.
  - *hopla_phasing.py*: to run phasing for all chromosomes, with the help
                        of Hopla to instanciate a *phase.py* call for
                        each chromosome (or X region).


=================================== ================================================================================= ===========
 *phase.py*
---------------------------------------------------------------------------------------------------------------------------------
            Parameter                                            Description                                           Default
=================================== ================================================================================= ===========
``-i, --bfile``                     Path without extension to the preprocessed dataset in Plink bed/bim/fam format.   *required*
``-c, --chromosome <name>``         Name of the chromosme or X region {1..22, X_PAR1, X_nonPAR, X_PAR2}.              *required*
``-p, --ref-hap <path>``            Path to reference panel haplotypes in Impute2 format.                             *required*
``-l, --ref-legend <path>``         Path to reference panel legend in Impute2 format.                                 *required*
``-s, --ref-sample <path>``         Path to reference panel sample in Impute2 format.                                 *required*
``-r, --recombination-map <path>``  Path to the recombination map of the chromosome in Impute2 format.                *required*
``-o, --outdir <path>``             Path to directory where to output.                                                *required*
``-b, --build <version>``           Genomic build version of the preprocessed/reference panel data (e.g. hg19).       *required*
``-n, --nb-cpus <int>``             Number of CPUs to use in the phasing step (Shapeit2 ``--thread`` argument).       *required*
``-z, --Ne <int>``                  *Shapeit2* ``--effective-size`` parameter. See also Impute2 ``-Ne`` parameter.    *20000*
``-j, --plink-exe <exe>``           Path to the Plink executable or alias if it's in $PATH.                           *'plink'*
``-y, --shapeit2-exe <exe>``        Path to the Shapeit2 executable or alias if it's in $PATH.                        *'shapeit'*
``-Z, --no-verbose``                Do not log.                                                                       *False*
=================================== ================================================================================= ===========

|

=================================== ================================================================================ ===========
 *hopla_phasing.py*
--------------------------------------------------------------------------------------------------------------------------------
            Parameter                                           Description                                            Default
=================================== ================================================================================ ===========
``-i, --bfile``                     Path without extension to the preprocessed dataset in Plink bed/bim/fam format.  *required*
``-p, --ref-panel <path>``           | Path to the table file describing the reference panel.                        *required*
                                     | This file is automatically created if the reference panel was created using
                                     | the provided script.
                                     | It relates a chromosome name to 4 reference files in Impute2 format:
                                     | reference haplotypes/legend/sample and recombination map.
                                     | The paths can be absolute or relative to the table file. File structure:
                                     | one row per chromosome, tab-separated fields
                                     | ``<CHROM>\t<REF HAP>\t<REF LEGEND>\t<REF SAMPLE>\t<RECOMBINATION MAP>``
                                     | Accepted chromosome names: {1..22, X_PAR1, X_nonPAR, X_PAR2}
``-o, --outdir <path>``             Path to directory where to output.                                               *required*
``-b, --build <version>``           Genomic build version of the preprocessed/reference panel data (e.g. hg19).      *required*
``-n, --nb-processes <int>``        Number of simultaneous *phase.py* call.                                          *required*
``-c, --nb-cpus-per-process <int>`` Value of ``--nb-cpus`` for each *phase.py* call.                                 *required*
``-z, --Ne <int>``                  Shapeit2 ``--effective-size`` parameter. See also Impute2 ``-Ne`` parameter.     *20000*
``-j, --plink-exe <exe>``           Path to the Plink executable or alias if it's in $PATH.                          *'plink'*
``-y, --shapeit2-exe <exe>``        Path to the Shapeit2 executable or alias if it's in $PATH.                       *'shapeit'*
=================================== ================================================================================ ===========


Strand alignment
~~~~~~~~~~~~~~~~
Shapeit2 provides an option ``-check`` that detects strand alignment
problems between the dataset to impute and the reference panel.
We use it to align the dataset, using the following recipe:
  - check strand alignment with Shapeit2
  - run Plink ``--flip`` to flip all variants that were detected as not aligned
  - re-check strand alignment with Shapeit2
  - remove variants that are still not aligned: they are considered as unalignable

