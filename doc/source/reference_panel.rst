
Reference panel
===============

A reference panel is required for both phasing and imputation. You can use
any reference panel as long as it is in the Shapeit2/Impute2 format.

The 1000 Genomes phase 3 reference panel is publicly available and can be
downloaded at `<https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference>`_

Alternatively, the package provides a script to create a reference panel
to impute datasets of european samples:

============================= ============================================================================ ============
*create_ref_panel_of_european_biallelic_variants.py*
-----------------------------------------------------------------------------------------------------------------------
        Parameter                                           Description                                       Default
============================= ============================================================================ ============
``-o, --outdir <path>``       Path to directory where to download, filter and convert the reference panel. *required*
``-n, --nb-processes <int>``  Number of processes to run in parallel (number of CPUs to use).              *required*
                              Each process filters and converts one chromosome of the reference panel.
``-e, --bcftools-exe <exe>``  Path to the Bcftools executable or alias if it's in $PATH.                   *'bcftools'*
============================= ============================================================================ ============

Because most of the variants in the 1000 Genomes phase 3 are monomorphic
in the european population and because imputation is computationally
expensive, we define a smaller reference panel for imputing european
samples.

Script steps:
  - download the 1000 Genomes phase 3 reference panel in VCF format
  - download the recombination maps in Impute2 format,
  - filter all the VCF files by keeping only variants that:
        - are bi-allelic,
        - of type SNP or INDEL
        - have a minor allele count >= 2 in the EUR super-population
          (~500 samples).
  - convert to Shapeit2/Impute2 format,

Notes:
  * only variants are filtered, all subjects from all super-populations
    are kept in the reference panel. Keeping all ancestries improves
    imputation of rare variants.
  * the X chromosome is splitted in 3 regions: X_PAR1, X_nonPAR and X_PAR2.
  * Downloading and filtering the ~17GB takes hours.

This script was inspired by the ENIGMA protocol where they used a smaller
custom reference panel from the 1000 Genomes phase 1 data:
http://enigma.ini.usc.edu/wp-content/uploads/2012/07/ENIGMA2_1KGP_cookbook_v3.pdf
