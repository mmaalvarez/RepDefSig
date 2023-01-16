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

// the same median_scores is passed to 2 diff processes
median_scores.into{median_scores_for_next_process; median_scores_for_last_process}



// in parallel per chromosome (no chrY, as RepliSeq does not have it), load genomic coordinates of the DNA repair marks and chromatin features that are specified in input_lists/

chromosomes = Channel.from( ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] )

process load_feature_maps {

    time = 5.hour
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }
    //memory = { (chromosomes.val>10) ? (params.memory_process2 + 5*(task.attempt-1)).GB : (params.memory_process2 - 20 + 5*(task.attempt-1)).GB }

    input:
    path 'dnarep_marks' from params.dnarep_marks
    path 'chromatin_features' from params.chromatin_features
    val 'chromosome' from chromosomes
    path median_scores_collected from median_scores_for_next_process.collect() // from previous process

    output:
    file 'map_features_chr*.tsv' into map_features
    file 'map_features_binarized_chr*.tsv' into map_features_binarized

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${median_scores_collected}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${median_scores_collected}
    fi
    """
}

// rowbind all chromosomes' raw score map features, will be read by final process's R script (hardcoded)
map_features
    .collectFile(name: "res/map_features.tsv", keepHeader: true)



// calculate offset (log(#trinuc+1)) for each of the 32 N(C|T)N types (e.g. ACT) that exist in each RT×dnarepmarks combination, and could therefore be any A(C>D)T SNV

process offset {

    time = 1.hour
    memory = { (params.memory_process3 + 5*(task.attempt-1)).GB }

    input:
    path map_features_binarized_collected from map_features_binarized.collect() // from previous process

    output:
    file 'offset.tsv' into offset_table

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/3_offset.R ${map_features_binarized_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/3_offset.R ${map_features_binarized_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/3_offset.R ${map_features_binarized_collected}
    fi
    """
}



// process in parallel each sample specified in input_lists/sample_ids.tsv

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header:false)
    .set{ sample_name }

process load_sample_somatic_muts_overlap_feature_maps_run_regressions {

    time = 6.hour
    memory = { (params.memory_process4 + 5*(task.attempt-1)).GB }

    input:
    val sample from sample_name
    val somatic_data from params.somatic_data
    val metadata from params.metadata
    path dnarep_marks from params.dnarep_marks
    path chromatin_features from params.chromatin_features 
    path offset from offset_table // from previous process
    path median_scores_collected from median_scores_for_last_process.collect() // from 1st process
    // res/map_features.tsv from 2nd process read directly in R 

    output:
    file 'results_sample.tsv' into results

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${median_scores_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${median_scores_collected}
    fi
    """
}

// rowbind regression results of all samples
results
    .collectFile(name: 'res/results.tsv', keepHeader: true)
    .println { "Finished! Combined results for all samples saved in res/results.tsv" }
