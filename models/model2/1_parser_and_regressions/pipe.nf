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
    path 'low_mappability_regions' from params.low_mappability_regions

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
            Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path} ${low_mappability_regions}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path} ${low_mappability_regions}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/1_list_medians_scores.R ${name} ${path} ${low_mappability_regions}
    fi
    """
}

// the same median_scores is passed to 3 diff processes
median_scores.into{median_scores_for_load_feature_maps; median_scores_for_load_sample_somatic_muts_overlap_feature_maps_run_regressions; median_scores_for_sim_pos_con}



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
    path 'low_mappability_regions' from params.low_mappability_regions
    path median_scores_collected from median_scores_for_load_feature_maps.collect() // from previous process

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
            Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${low_mappability_regions} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${low_mappability_regions} ${median_scores_collected}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_load_feature_maps.R ${dnarep_marks} ${chromatin_features} ${chromosome} ${low_mappability_regions} ${median_scores_collected}
    fi
    """
}

// rowbind all chromosomes' raw score map features, will be read by 4th process's R script (hardcoded)
map_features
    .collectFile(name: "res/map_features.tsv", keepHeader: true)



// calculate offset (log(#trinuc+1)) for each of the 32 N(C|T)N types (e.g. ACT) that exist in each RT×dnarepmarks combination, and could therefore be any A(C>D)T SNV

process offset {

    time = 1.hour
    memory = { (params.memory_process3 + 5*(task.attempt-1)).GB }

    input:
    path 'low_mappability_regions' from params.low_mappability_regions
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
            Rscript $PWD/bin/3_offset.R ${low_mappability_regions} ${map_features_binarized_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/3_offset.R ${low_mappability_regions} ${map_features_binarized_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/3_offset.R ${low_mappability_regions} ${map_features_binarized_collected}
    fi
    """
}

// the same offset_table is passed to 2 diff processes
offset_table.into{offset_table_for_load_sample_somatic_muts_overlap_feature_maps_run_regressions; offset_table_for_sim_pos_con}



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
    path offset from offset_table_for_load_sample_somatic_muts_overlap_feature_maps_run_regressions // from previous process
    path 'low_mappability_regions' from params.low_mappability_regions
    path median_scores_collected from median_scores_for_load_sample_somatic_muts_overlap_feature_maps_run_regressions.collect() // from 1st process
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
            Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${median_scores_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${median_scores_collected}
    fi
    """
}

// rowbind regression results of all (real) samples
results
    .collectFile(name: 'res/results.tsv', keepHeader: true)
    .println { "Regression results for all real samples saved in res/results.tsv" }



// simulate positive controls and run regressions for them

mutation_foldinc = Channel.from( params.mutation_foldinc ) // in parallel, fold-values by which increase mutation burden in each DNA repair mark's "high abundance" bins

dnarep_marks_simulate = Channel.from( EXTRACT_NAMES(params.dnarep_marks) ) // also in parallel, which dna repair mark has mutations increased

process sim_pos_con {
    
    time = 6.hour
    memory = { (params.memory_process5 + 5*(task.attempt-1)).GB }

    input:
    val somatic_data from params.somatic_data
    val metadata from params.metadata
    path dnarep_marks from params.dnarep_marks
    path chromatin_features from params.chromatin_features 
    path offset from offset_table_for_sim_pos_con // from 3rd process
    path 'low_mappability_regions' from params.low_mappability_regions
    val 'mutation_foldinc' from mutation_foldinc
    val 'dnarep_marks_simulate' from dnarep_marks_simulate
    path median_scores_collected from median_scores_for_sim_pos_con.collect() // from 1st process

    output:
    file 'simulated_positive_control.tsv' into sim_pos_con

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/5_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${mutation_foldinc} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/5_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${mutation_foldinc} ${median_scores_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/5_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} $PWD ${low_mappability_regions} ${mutation_foldinc} ${dnarep_marks_simulate} ${median_scores_collected}
    fi
    """
}

// rowbind regression results of all simulated samples
sim_pos_con
    .collectFile(name: 'res/simulated_positive_controls.tsv', keepHeader: true)
    .println { "Regression results for simulated positive control samples saved in res/simulated_positive_controls.tsv" }


println "Finished!"
