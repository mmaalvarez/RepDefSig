#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R SETD2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/DNA_repair/MMR/MSH6_SETD2_H3K36me3_guominli/2_fold_enrichment_vs_input/bedgraphs_fold_enrich/SRR13314134_sorted_FE.bdg
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R SETD2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/DNA_repair/MMR/MSH6_SETD2_H3K36me3_guominli/2_fold_enrichment_vs_input/bedgraphs_fold_enrich/SRR13314134_sorted_FE.bdg
    fi
else
    # no conda installed
    Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R SETD2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/DNA_repair/MMR/MSH6_SETD2_H3K36me3_guominli/2_fold_enrichment_vs_input/bedgraphs_fold_enrich/SRR13314134_sorted_FE.bdg
fi
