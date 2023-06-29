#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --dnarep_marks $PWD/input_lists/dnarep_marks.csv \
												--chromatin_features $PWD/input_lists/chromatin_features.csv \
												--good_mappability_regions /g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed \
												--utils /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R \
												--memory_process1 13 \
												--memory_process2 30 \
												--memory_process3 80 \
												--memory_process4 15 \
												-resume #\ -with-tower
