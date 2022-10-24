#!/bin/bash

export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --dnarep_marks $PWD/input_lists/dnarep_marks.tsv \
												--chromatin_features $PWD/input_lists/chromatin_features.tsv \
												--sample_ids $PWD/input_lists/sample_ids.tsv \
												--somatic_data /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_ \
												--metadata /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv \
												--memory_process1 2 \
												--memory_process2 2 \
												-with-tower \
												-resume #-bg
