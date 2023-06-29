#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/3_binarize_scores.R chromatin_features map_features_single_chr median_score_DHS.tsv median_score_RnaSeq.tsv median_score_exons.tsv median_score_H3K36me3.tsv median_score_RepliSeq.tsv median_score_SETD2.tsv median_score_MSH6.tsv
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/3_binarize_scores.R chromatin_features map_features_single_chr median_score_DHS.tsv median_score_RnaSeq.tsv median_score_exons.tsv median_score_H3K36me3.tsv median_score_RepliSeq.tsv median_score_SETD2.tsv median_score_MSH6.tsv
    fi
else
    # no conda installed
    Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/3_binarize_scores.R chromatin_features map_features_single_chr median_score_DHS.tsv median_score_RnaSeq.tsv median_score_exons.tsv median_score_H3K36me3.tsv median_score_RepliSeq.tsv median_score_SETD2.tsv median_score_MSH6.tsv
fi
