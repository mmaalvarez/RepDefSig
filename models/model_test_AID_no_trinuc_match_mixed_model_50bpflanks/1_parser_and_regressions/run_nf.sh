#!/bin/bash

# export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --AID_regions $PWD/input_lists/AID_regions.csv \
												--chromatin_features $PWD/input_lists/chromatin_features.csv \
												--sample_ids $PWD/input_lists/sample_ids.csv \
												--somatic_data /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_ \
												--low_mappability_regions /g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed \
												--memory_process2 3 \
												--memory_process3 2 \
												--memory_process4 3 \
												-resume #\ -with-tower
