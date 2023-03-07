#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --dnarep_marks $PWD/input_lists/dnarep_marks.csv \
												--chromatin_features $PWD/input_lists/chromatin_features.csv \
												--sample_ids $PWD/input_lists/sample_ids.csv \
												--somatic_data /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_ \
												--metadata /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv \
												--mutation_foldinc 2,4,8,16 \
												--good_mappability_regions /g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed \
												--utils /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R \
												--memory_process1 10 \
												--memory_process2 15 \
												--memory_process3 100 \
												--memory_process4 8 \
												--memory_process5.1 15 \
												--memory_process5.2 15 \
												--memory_process6.1 15 \
												--memory_process6.2 15 \
												-resume #\ -with-tower
