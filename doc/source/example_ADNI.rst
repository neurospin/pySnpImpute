
Imputation from scratch: example with the ADNI dataset
======================================================

Example of imputation of ADNI genotypes:
  - using imputed sex (SNP sex) instead of reported sex
  - assuming one multicore computer with 20 available CPUs

.. code-block:: bash

  ## Setting up environment/data variables

  # Directory where to download and format the reference panel
  ref_panel_root_dir=/volatile/imputation/references/1000GP_phase3

  # Path without extension to the bed/bim/fam to impute, here ADNI data
  raw_data=/volatile/imputation/ADNI/raw_data/ADNI_cluster_01_forward_757LONI

  # Build versions
  raw_data_build=hg18
  ref_panel_build=hg19

  # Effective size of population, see also Impute2 Ne parameter
  Ne=20000

  # Directory where to output
  outdir=/volatile/imputation/ADNI

  # Assuming we are running on one multicore computer
  nb_cpus=20

  # Assuming you have installed the package, get /scripts/ directory
  scripts_dir=`python -c "import pysnpimpute; print pysnpimpute.__path__[0]"`/scripts

  # --------------------------------------------------------------------
  # Run

  # 1 - Create the reference panel
  python ${scripts_dir}/create_ref_panel_of_european_biallelic_variants.py \
        --outdir ${ref_panel_root_dir} \
        --nb-processes ${nb_cpus}

  # Generated file describing the reference panel
  ref_panel=${ref_panel_root_dir}/Impute2_european_variants/ref_panel.txt

  # 2 - Preprocess
  preproc_outdir=${outdir}/preprocessed_data
  python ${scripts_dir}/preprocess.py \
        --bfile ${raw_data} \
        --outdir ${preproc_outdir} \
        --lift ${raw_data_build} ${ref_panel_build} \
        --use-imputed-sex

  # 3 - Phase
  phasing_outdir=${outdir}/phasing
  python ${scripts_dir}/hopla_phasing.py \
        --bfile ${preproc_outdir}/ADNI_cluster_01_forward_757LONI.hg19.imputedsex.qc \
        --ref-panel ${ref_panel} \
        --outdir ${phasing_outdir} \
        --build ${ref_panel_build} \
        --Ne ${Ne} \
        --nb-processes 1 \
        --nb-cpus-per-process ${nb_cpus}

  # Generated file describing the phased dataset
  phased_dataset=${phasing_outdir}/phased_dataset.txt

  # 4 - Impute
  python ${scripts_dir}/hopla_imputation.py \
        --phased-dataset ${phased_dataset} \
        --ref-panel ${ref_panel} \
        --outdir ${outdir}/imputation \
        --build ${ref_panel_build} \
        --Ne ${Ne} \
        --nb-processes ${nb_cpus}
