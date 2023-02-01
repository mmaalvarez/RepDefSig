#!/usr/bin/env nextflow


// in parallel per chromosome (no chrY, as RepliSeq does not have it), load genomic coordinates of the AID regions and RepTime, specified in input_lists/

chromosomes = Channel.from( ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] )

process load_feature_maps {

    time = 1.hour
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    path 'AID_regions' from params.AID_regions
    path 'chromatin_features' from params.chromatin_features
    val 'chromosome' from chromosomes
    path 'low_mappability_regions' from params.low_mappability_regions

    output:
    file 'map_features_chr*.tsv' into map_features
    file 'map_features_binarized_chr*.tsv' into map_features_binarized

    """
    #!/usr/bin/env bash

    Rscript $PWD/bin/2_load_feature_maps.R ${AID_regions} ${chromatin_features} ${chromosome} ${low_mappability_regions}
    """
}

// rowbind all chromosomes' raw score map features, will be read by last process's R script (hardcoded)
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

    Rscript $PWD/bin/3_offset.R ${low_mappability_regions} ${map_features_binarized_collected}
    """
}



// process in parallel each sample specified in input_lists/sample_ids.tsv

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header:false)
    .set{ sample_name }

process load_sample_somatic_muts_overlap_feature_maps_run_regressions {

    time = 1.hour
    memory = { (params.memory_process4 + 5*(task.attempt-1)).GB }

    input:
    val sample from sample_name
    val somatic_data from params.somatic_data
    path AID_regions from params.AID_regions
    path chromatin_features from params.chromatin_features 
    path offset from offset_table // from previous process
    path 'low_mappability_regions' from params.low_mappability_regions
    // res/map_features.tsv from 1st process read directly in R 

    output:
    file 'results_sample.tsv' into results

    """
    #!/usr/bin/env bash

    Rscript $PWD/bin/4_load_sample_somatic_muts_overlap_feature_maps_run_regressions.R ${sample} ${somatic_data} ${AID_regions} ${chromatin_features} ${offset} $PWD ${low_mappability_regions}
    """
}

// rowbind regression results of all samples
results
    .collectFile(name: 'res/results.tsv', keepHeader: true)
    .println { "Finished! Regression results for all real samples saved in res/results.tsv" }
