#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/4_offset.R map_features_binarized_chr21.tsv map_features_binarized_chr22.tsv map_features_binarized_chr19.tsv map_features_binarized_chr20.tsv map_features_binarized_chr18.tsv map_features_binarized_chr14.tsv map_features_binarized_chr16.tsv map_features_binarized_chr13.tsv map_features_binarized_chr9.tsv map_features_binarized_chrX.tsv map_features_binarized_chr12.tsv map_features_binarized_chr17.tsv map_features_binarized_chr10.tsv map_features_binarized_chr15.tsv map_features_binarized_chr7.tsv map_features_binarized_chr4.tsv map_features_binarized_chr11.tsv map_features_binarized_chr8.tsv map_features_binarized_chr6.tsv map_features_binarized_chr3.tsv map_features_binarized_chr2.tsv map_features_binarized_chr5.tsv map_features_binarized_chr1.tsv
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/4_offset.R map_features_binarized_chr21.tsv map_features_binarized_chr22.tsv map_features_binarized_chr19.tsv map_features_binarized_chr20.tsv map_features_binarized_chr18.tsv map_features_binarized_chr14.tsv map_features_binarized_chr16.tsv map_features_binarized_chr13.tsv map_features_binarized_chr9.tsv map_features_binarized_chrX.tsv map_features_binarized_chr12.tsv map_features_binarized_chr17.tsv map_features_binarized_chr10.tsv map_features_binarized_chr15.tsv map_features_binarized_chr7.tsv map_features_binarized_chr4.tsv map_features_binarized_chr11.tsv map_features_binarized_chr8.tsv map_features_binarized_chr6.tsv map_features_binarized_chr3.tsv map_features_binarized_chr2.tsv map_features_binarized_chr5.tsv map_features_binarized_chr1.tsv
    fi
else
    # no conda
    Rscript /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_MMR_poster/2_signature_extraction_and_plotting/1_QC/trinuc_proportions_unmatched_feature_tracks/bin/4_offset.R map_features_binarized_chr21.tsv map_features_binarized_chr22.tsv map_features_binarized_chr19.tsv map_features_binarized_chr20.tsv map_features_binarized_chr18.tsv map_features_binarized_chr14.tsv map_features_binarized_chr16.tsv map_features_binarized_chr13.tsv map_features_binarized_chr9.tsv map_features_binarized_chrX.tsv map_features_binarized_chr12.tsv map_features_binarized_chr17.tsv map_features_binarized_chr10.tsv map_features_binarized_chr15.tsv map_features_binarized_chr7.tsv map_features_binarized_chr4.tsv map_features_binarized_chr11.tsv map_features_binarized_chr8.tsv map_features_binarized_chr6.tsv map_features_binarized_chr3.tsv map_features_binarized_chr2.tsv map_features_binarized_chr5.tsv map_features_binarized_chr1.tsv
fi
