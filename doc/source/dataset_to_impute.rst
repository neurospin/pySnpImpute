
Dataset to impute
=================

Requirements regarding the dataset to impute:
  - has to be provided as a Plink binary dataset (bed/bim/fam)
  - made of unrelated samples (subjects)
  - the X chromosome variants should be splitted according the to Plink
    convention, i.e. the variants of the PseudoAutosomal regions of the
    X chromosome should be associated to chromosome 'XY' (or 25) and not
    to 'X' (or 23). This can be handled with the ``--split-x`` option
    of Plink.

.. The dataset to impute has to be provided as a Plink binary dataset (bed/bim/fam),
   made of unrelated samples (subjects) and the X chromosome variants should
   be splitted according the to Plink convention, i.e. the variants of the
   PseudoAutosomal regions of the X chromosome should be associated to
   chromosome 'XY' (or 25) and not to 'X' (or 23). This can be handled with
   the --split-x option of Plink.
