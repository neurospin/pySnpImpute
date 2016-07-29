
Imputation
==========

Chromosomes are not imputed as a whole but divided in chunks that are
imputed separately. Like the phasing step, the package provides 2 scripts
to run imputation.

The first script runs imputations on a single chunk and can be used to
impute a small region (typically < 5Mb).
The second scripts runs a complete imputation. It divides the chromosomes
into chunks and runs the first script for each one of them.

=================================== =============================================================================== ===========
 *impute.py*
-------------------------------------------------------------------------------------------------------------------------------
            Parameter                                            Description                                         Default
=================================== =============================================================================== ===========
``-c, --chromosome <name>``         Name of the chromosme or X region {1..22, X_PAR1, X_nonPAR, X_PAR2}.            *required*
``-f, --from-bp``                   Basepair position where to start imputation.                                    *required*
``-t, --to-bp``                     Basepair position where to stop imputation.                                     *required*
``-i, --hap``                       Path to haplotypes to impute in Impute2 format.                                 *required*
``-s, --sample``                    Path to sample file associated to haplotypes to impute in Impute2 format.       *required*
``-p, --ref-hap <path>``            Path to reference panel haplotypes in Impute2 format.                           *required*
``-l, --ref-legend <path>``         Path to reference panel legend in Impute2 format.                               *required*
``-r, --recombination-map <path>``  Path to the recombination map of the chromosome in Impute2 format.              *required*
``-o, --outdir <path>``             Path to directory where to output.                                              *required*
``-b, --build <version>``           Genomic build version of the <bfile>/ref panel data (e.g. hg19).                *required*
``-n, --nb-cpus <int>``             Number of CPUs to use in the phasing step (Shapeit2 ``--thread`` argument).     *required*
``-z, --Ne <int>``                  Impute2 ``-Ne`` parameter.                                                      *20000*
``-u, --buffer <int>``              Impute2 ``-buffer`` parameter. Length of buffer region in Kb (NOT IN BASEPAIR)  *250*
                                    to include on each side of the analysis interval.                     
``-g, --allow-large-regions``       By default Impute2 does not allow imputation on regions of size > 7 Mb.         *False*
                                    To force imputation on bigger regions, set this option. Note that the buffers,
                                    on each side of the imputation region, have be taken into account.
``-e, --impute2-exe <exe>``         Path to the Impute2 executable or alias if it's in $PATH.                       *'impute2'*
``-z, --no-verbose``                Do not log.                                                                     *False*
=================================== =============================================================================== ===========

|

=================================== =============================================================================== ===========
 *hopla_imputation.py*
-------------------------------------------------------------------------------------------------------------------------------
            Parameter                                            Description                                         Default
=================================== =============================================================================== ===========
``-i, --phased_dataset``            | Path to a table file listing the paths to the phased haplotypes to impute     *required*
                                    | (in Impute2 format) and the related chromosome names. 
                                    | This file is automatically created if the phasing was done with
                                    | *hopla_phasing.py*. 
                                    | The paths can be absolute or relative to the table file.
                                    | File structure: One row per input file, tab-separated fields:
                                    | ``<CHROM>\t<PHASED HAP>\t<SAMPLE>``
                                    | Accepted chromosome names: {1..22, X_PAR1, X_nonPAR, X_PAR2}
``-p, --ref-panel <path>``          | Path to the table file describing the reference panel.                        *required*
                                    | This file is automatically created if the reference panel was created using
                                    | the provided script.
                                    | It relates a chromosome name to 4 reference files in Impute2 format:
                                    | reference haplotypes/legend/sample and recombination map.
                                    | The paths can be absolute or relative to the table file. File structure:
                                    | one row per chromosome, tab-separated fields
                                    | ``<CHROM>\t<REF HAP>\t<REF LEGEND>\t<REF SAMPLE>\t<RECOMBINATION MAP>``
                                    | Accepted chromosome names: {1..22, X_PAR1, X_nonPAR, X_PAR2}
``-o, --outdir <path>``             Path to directory where to output.                                              *required*
``-b, --build <version>``           Genomic build version of the to-impute/ref haplotypes (e.g. hg19).              *required*
``-n, --nb-processes <int>``        Each process runs the imputation of one chunk with *impute.py*                  *required*
``-z, --Ne <int>``                  Impute2 ``-Ne`` parameter.                                                      *20000*
``-k, --chunksize-kb <int>``        The chromosome is imputed by chunks of size ~<chunksize_kb> Kb.                 *5000*
                                    By default 5000 (i.e. 5Mb).
``-N <int>``                        The minimum number of known variants per chunk. A chunk will cover an interval  *200*
                                    of length <chunksize> Kb except if there is less than <N> variants in the
                                    interval. In that case the chunk is enlarged to ensure that there is at least
                                    N variants. By default 200.
``-u, --buffer <int>``              Impute2 ``-buffer`` parameter. Length of buffer region in Kb (NOT IN BASEPAIR)  *250*
                                    to include on each side of the analysis interval.

``-e, --impute2-exe <exe>``         Path to the Impute2 executable or alias if it's in $PATH.                       *'impute2'*
``-z, --no-verbose``                Do not log.                                                                     *False*
=================================== =============================================================================== ===========


Chunking of chromosomes
~~~~~~~~~~~~~~~~~~~~~~~

For each chromosome *hopla_imputation.py* computes intervals that define
chromosome chunks to impute separately.

It does so by following a few rules:
    - the chunks should not overlap with each other
    - the chunks should not overlap with the centromere region of the chromosome.
      Each chromosome arm is chunked separately.
    - each chunk should have at least ``--N <N>`` known haplotypes. By default 200.
    - each chunk should have a minimum size of ``--chunksize <chunksize_kb>`` Kb.
      The chunk is made bigger if there not enough known haplotypes in the interval.

If a variant cannot be placed inside a chunk that respects these constraints,
then it is not imputable or, for known variants, not usable for imputation.
