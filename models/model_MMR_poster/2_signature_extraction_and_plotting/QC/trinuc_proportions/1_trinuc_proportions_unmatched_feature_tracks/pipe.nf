#!/usr/bin/env nextflow



// calculate median score of each dna rep mark in parallel

Channel
    .fromPath(params.dnarep_marks)
    .splitCsv(header:true)
    .map{ row-> tuple(row.name, row.path) }
    .set{ dnarep_mark_paths }

process list_median_scores {

    time = 2.hour
    memory = { (params.memory_process1 + 5*(task.attempt-1)).GB }

    input:
    set name, path from dnarep_mark_paths

    output:
    file 'median_score_*.tsv' into median_scores

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path}
    fi
    """
}



// in parallel per chromosome (no chrY, as RepliSeq does not have it), load genomic coordinates of the DNA repair marks and chromatin features that are specified in input_lists/

chromosomes = Channel.from( ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] )

process load_feature_maps {

    time = 2.hour
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    path 'utils' from params.utils
    path 'dnarep_marks' from params.dnarep_marks
    path 'chromatin_features' from params.chromatin_features
    val 'chromosome' from chromosomes

    output:
    file 'map_features_chr*.tsv' into map_features

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/2_load_feature_maps.R ${utils} ${dnarep_marks} ${chromatin_features} ${chromosome}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_load_feature_maps.R ${utils} ${dnarep_marks} ${chromatin_features} ${chromosome}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_load_feature_maps.R ${utils} ${dnarep_marks} ${chromatin_features} ${chromosome}
    fi
    """
}


process binarize_scores {

    time = 9.hour
    memory = { (params.memory_process3 + 40*(task.attempt-1)).GB }

    input:
    path 'chromatin_features' from params.chromatin_features
    path 'map_features_single_chr' from map_features // raw score map features PER CHROMOSOME (not collected)
    path median_scores_collected from median_scores.collect() // from first process

    output:
    file 'map_features_binarized_chr*.tsv' into map_features_binarized

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/3_binarize_scores.R ${chromatin_features} ${map_features_single_chr} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/3_binarize_scores.R ${chromatin_features} ${map_features_single_chr} ${median_scores_collected}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/3_binarize_scores.R ${chromatin_features} ${map_features_single_chr} ${median_scores_collected}
    fi
    """
}



// calculate offset (log(#trinuc+1)) for each of the 32 N(C|T)N types (e.g. ACT) that exist in each RT×dnarepmarks combination, and could therefore be any A(C>D)T SNV

process offset {

    publishDir "$PWD/res/", mode: 'copy'

    time = 1.hour
    memory = { (params.memory_process4 + 5*(task.attempt-1)).GB }

    input:
    path map_features_binarized_collected from map_features_binarized.collect() // from previous process

    output:
    file 'offset.tsv'
    file 'trinuc_dist_QC.jpg'

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/4_offset.R ${map_features_binarized_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/4_offset.R ${map_features_binarized_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/4_offset.R ${map_features_binarized_collected}
    fi
    """
}
