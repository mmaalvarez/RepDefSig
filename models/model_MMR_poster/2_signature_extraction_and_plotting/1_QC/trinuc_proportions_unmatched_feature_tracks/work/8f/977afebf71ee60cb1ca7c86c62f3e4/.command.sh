#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R RepliSeq /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/replication_time/fran_pooled_bins_replication_time/fran_legacy_parsed/RepliSeq_pooled8_merged_1-3_vs_4-6.bed
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R RepliSeq /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/replication_time/fran_pooled_bins_replication_time/fran_legacy_parsed/RepliSeq_pooled8_merged_1-3_vs_4-6.bed
    fi
else
    # no conda installed
    Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R RepliSeq /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/replication_time/fran_pooled_bins_replication_time/fran_legacy_parsed/RepliSeq_pooled8_merged_1-3_vs_4-6.bed
fi
