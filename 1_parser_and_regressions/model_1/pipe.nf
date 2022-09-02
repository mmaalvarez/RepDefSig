#!/usr/bin/env nextflow


// this loads genomic coordinates of the DNA repair marks and chromatin features that are specified in input_lists/

process load_feature_maps {

    time = 1.hour
    memory = { (params.memory_process1 + 5*(task.attempt-1)).GB }

    input:
    path 'dnarep_marks' from params.dnarep_marks
    path 'chromatin_features' from params.chromatin_features

    output:
    file 'ln_bpsum_chromatin_env_table.tsv' into ln_bpsum_chromatin_env_table
    file 'map_features.tsv' into map_features

    """
    #!/usr/bin/env bash

    conda activate R

    Rscript $PWD/bin/1_load_feature_maps.R ${dnarep_marks} ${chromatin_features}
    """
}


// process in parallel each sample specified in input_lists/sample_ids.tsv

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header:false)
    .set{ sample_name }

process load_sample_somatic_muts_overlap_feature_maps_run_regressions {

    time = 1.hour
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    // these are specified in the command parameters
    val sample from sample_name
    val somatic_data from params.somatic_data
    path metadata from params.metadata
    path dnarep_marks from params.dnarep_marks
    path chromatin_features from params.chromatin_features 

    // these are read directly in the R script -- it's the output from the previous process
    path 'ln_bpsum_chromatin_env_table.tsv' from ln_bpsum_chromatin_env_table
    path 'map_features.tsv' from map_features

    output:
    file 'results_sample.tsv' into results

    """
    #!/usr/bin/env bash

    conda activate R

    Rscript $PWD/bin/2_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${metadata} ${dnarep_marks} ${chromatin_features} 
    """
}


// combine regression results of all samples

results
    .collectFile(name: 'res/results.tsv', keepHeader: true)
    .println { "Finished! Combined results for all samples saved in res/results.tsv" }
