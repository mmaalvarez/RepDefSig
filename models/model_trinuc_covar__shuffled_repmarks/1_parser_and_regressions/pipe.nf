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

// the same median_scores is passed to several processes
median_scores.into{median_scores_for_binarize_scores; median_scores_for_load_sample_somatic_muts_overlap_feature_maps; median_scores_for_sim_pos_con}



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

// the same map_features is passed to several processes
map_features.into{map_features_for_binarize_scores; map_features_for_process5_1; map_features_for_process6_1}



process binarize_scores {

    time = 9.hour
    memory = { (params.memory_process3 + 20*(task.attempt-1)).GB }

    input:
    path 'chromatin_features' from params.chromatin_features
    path 'map_features_single_chr' from map_features_for_binarize_scores // raw score map features PER CHROMOSOME (not collected)
    path median_scores_collected from median_scores_for_binarize_scores.collect() // from first process

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

    time = 1.hour
    memory = { (params.memory_process4 + 5*(task.attempt-1)).GB }

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

// the same offset_table is passed to 2 diff processes
offset_table.into{offset_table_for_process5_2; offset_table_for_process6_1}




// process in parallel each sample specified in input_lists/sample_ids.tsv
// also in parallel per chromosome

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header:false)
    .set{ sample_name }

sample_name.into{sample_name_for_5_1; sample_name_for_5_2}

process load_sample_somatic_muts_overlap_feature_maps {

    time = 2.hour
    memory = { (params.memory_process5_1 + 4*(task.attempt-1)).GB }

    input:
    set sample, map_features_single_chr from sample_name_for_5_1.combine(map_features_for_process5_1) // nest raw score map features PER CHROMOSOME (not collected) within sample_name
    val somatic_data from params.somatic_data
    path dnarep_marks from params.dnarep_marks
    path chromatin_features from params.chromatin_features 
    path 'good_mappability_regions' from params.good_mappability_regions
    path median_scores_collected from median_scores_for_load_sample_somatic_muts_overlap_feature_maps.collect() // from 1st process

    output:
    file 'ready_for_regression_*_chr*.tsv' into ready_for_regression_single_chromosome

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/5.1_load_sample_somatic_muts_overlap_feature_maps.R ${sample} ${somatic_data} ${dnarep_marks} ${chromatin_features} ${map_features_single_chr} ${good_mappability_regions} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/5.1_load_sample_somatic_muts_overlap_feature_maps.R ${sample} ${somatic_data} ${dnarep_marks} ${chromatin_features} ${map_features_single_chr} ${good_mappability_regions} ${median_scores_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/5.1_load_sample_somatic_muts_overlap_feature_maps.R ${sample} ${somatic_data} ${dnarep_marks} ${chromatin_features} ${map_features_single_chr} ${good_mappability_regions} ${median_scores_collected}
    fi
    """
}



// now concatenate chromosomes and run regressions
process real_data_concat_chr_add_offset_run_regression {

    //queue = 'normal_prio_long'
    time = 8.hour
    memory = { (params.memory_process5_2 + 4*(task.attempt-1)).GB }

    input:
    path 'utils' from params.utils
    val sample from sample_name_for_5_2
    val metadata from params.metadata
    path dnarep_marks from params.dnarep_marks
    path offset from offset_table_for_process5_2 // from 4th process
    val trinuc_mode from params.trinuc_mode
    path ready_for_regression from ready_for_regression_single_chromosome.collect() // rowbind chromosomes (and all samples, each sample will be filtered in R script) from previous process

    output:
    file 'results_real_sample.tsv' into results_real_sample

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/5.2_concat_chr_add_offset_run_regression.R ${utils} ${sample} ${metadata} ${dnarep_marks} ${offset} ${trinuc_mode} ${ready_for_regression} 
        else
            # no conda environment named "R"
            Rscript $PWD/bin/5.2_concat_chr_add_offset_run_regression.R ${utils} ${sample} ${metadata} ${dnarep_marks} ${offset} ${trinuc_mode} ${ready_for_regression} 
        fi
    else
        # no conda
        Rscript $PWD/bin/5.2_concat_chr_add_offset_run_regression.R ${utils} ${sample} ${metadata} ${dnarep_marks} ${offset} ${trinuc_mode} ${ready_for_regression} 
    fi
    """
}

// rowbind regression results of all (real) samples
results_real_sample
    .collectFile(name: 'res/results_real_samples.tsv', keepHeader: true)
    .println { "Regression results for all real samples saved in res/results_real_samples.tsv" }





// simulate positive controls and run regressions for them
// also in parallel per chromosome

// in parallel, which dna repair mark has mutations increased
Channel
    .fromPath(params.dnarep_marks)
    .splitCsv(header:true)
    .map{ row-> tuple(row.name, row.path) }
    .set{ dnarep_marks_simulate }

// also in parallel, fold-values by which increase mutation burden in each DNA repair mark's "high abundance" bins
mutation_foldincs = Channel.from(params.mutation_foldinc.tokenize(','))

process sim_pos_con {
    
    time = 2.hour
    memory = { (params.memory_process6_1 + 5*(task.attempt-1)).GB }

    input:
    val somatic_data from params.somatic_data
    val metadata from params.metadata
    path dnarep_marks from params.dnarep_marks
    path chromatin_features from params.chromatin_features 
    path offset from offset_table_for_process6_1 // from 4th process
    path 'good_mappability_regions' from params.good_mappability_regions
    set name, path, mutation_foldinc, map_features_single_chr from dnarep_marks_simulate.combine(mutation_foldincs).combine(map_features_for_process6_1) // nest raw score map features PER CHROMOSOME (not collected) within mutation_foldincs within dnarep_marks_simulate
    path median_scores_collected from median_scores_for_sim_pos_con.collect() // from 1st process

    output:
    file 'ready_for_regression_sim_pos_con_chr*_mutfoldinc*_*.tsv' into ready_for_regression_single_chromosome_sim_pos_con

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/6.1_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} ${map_features_single_chr} ${good_mappability_regions} ${mutation_foldinc} ${name} ${median_scores_collected}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/6.1_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} ${map_features_single_chr} ${good_mappability_regions} ${mutation_foldinc} ${name} ${median_scores_collected}
        fi
    else
        # no conda
        Rscript $PWD/bin/6.1_simulate_pos_controls.R ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} ${offset} ${map_features_single_chr} ${good_mappability_regions} ${mutation_foldinc} ${name} ${median_scores_collected}
    fi
    """
}


// now concatenate chromosomes and run regressions

// AGAIN in parallel, which dna repair mark has mutations increased
Channel
    .fromPath(params.dnarep_marks)
    .splitCsv(header:true)
    .map{ row-> tuple(row.name, row.path) }
    .set{ dnarep_marks_simulate }

// AGAIN also in parallel, fold-values by which increase mutation burden in each DNA repair mark's "high abundance" bins
mutation_foldincs = Channel.from(params.mutation_foldinc.tokenize(','))

process sim_pos_con_concat_chr_run_regression {

    //queue = 'normal_prio_long'
    time = 8.hour
    memory = { (params.memory_process6_2 + 10*(task.attempt-1)).GB }

    input:
    path 'utils' from params.utils
    path dnarep_marks from params.dnarep_marks
    path offset from offset_table_for_process6_1 // from 4th process
    val trinuc_mode from params.trinuc_mode
    set name, path, mutation_foldinc from dnarep_marks_simulate.combine(mutation_foldincs) // nest mutation_foldincs within dnarep_marks_simulate
    path ready_for_regression_sim_pos_con from ready_for_regression_single_chromosome_sim_pos_con.collect() // combine chromosomes from previous process

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
            Rscript $PWD/bin/6.2_concat_chr_run_regression.R ${utils} ${dnarep_marks} ${name} ${mutation_foldinc} ${offset} ${trinuc_mode} ${ready_for_regression_sim_pos_con}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/6.2_concat_chr_run_regression.R ${utils} ${dnarep_marks} ${name} ${mutation_foldinc} ${offset} ${trinuc_mode} ${ready_for_regression_sim_pos_con}
        fi
    else
        # no conda
        Rscript $PWD/bin/6.2_concat_chr_run_regression.R ${utils} ${dnarep_marks} ${name} ${mutation_foldinc} ${offset} ${trinuc_mode} ${ready_for_regression_sim_pos_con}
    fi
    """
}

// rowbind regression results of all simulated samples
sim_pos_con
    .collectFile(name: 'res/simulated_positive_controls.tsv', keepHeader: true)
    .println { "Regression results for simulated positive control samples saved in res/simulated_positive_controls.tsv" }
