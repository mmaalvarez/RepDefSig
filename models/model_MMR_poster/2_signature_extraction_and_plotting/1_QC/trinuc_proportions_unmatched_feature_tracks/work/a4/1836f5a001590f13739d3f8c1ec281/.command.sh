#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R H3K36me2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/histone_methylation_and_DHS/H3K36me2/pooled_samples/pooled_all/collapse_ranges/H3K36me2_pooled_GSE118954-149670-175750.all_chr.collapsed_ranges.bed.gz
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R H3K36me2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/histone_methylation_and_DHS/H3K36me2/pooled_samples/pooled_all/collapse_ranges/H3K36me2_pooled_GSE118954-149670-175750.all_chr.collapsed_ranges.bed.gz
    fi
else
    # no conda installed
    Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/1_list_medians_scores.R H3K36me2 /g/strcombio/fsupek_cancer3/malvarez/chromatin_info/histone_methylation_and_DHS/H3K36me2/pooled_samples/pooled_all/collapse_ranges/H3K36me2_pooled_GSE118954-149670-175750.all_chr.collapsed_ranges.bed.gz
fi
