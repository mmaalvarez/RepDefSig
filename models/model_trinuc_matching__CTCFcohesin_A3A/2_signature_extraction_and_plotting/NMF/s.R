library(tidyverse)
library(data.table)
library(parallel)
library(msm) # rtnorm
library(NMF) # for NMF in first part (bootstrap) to determine optimal k signatures
library(RcppML) # for final NMF
library(lsa)
library(cluster)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("cosine", "lsa")
conflict_prefer("clusterExport", "parallel")
conflict_prefer("clusterEvalQ", "parallel")
conflict_prefer("parLapply", "parallel")
conflict_prefer("parLapplyLB", "parallel")
conflict_prefer("nmf", "NMF")


## name conversion tables 

PCAWG_names_conversion_table = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/pre_metadata/OLD_metadata.tsv") %>%
  rename("sample_id2" = "sample_id",
         "sample_id" = "icgc_donor_id") %>%
  filter(source %in% c("PCAWG")) %>% 
  mutate(icgc_donor_id = NA) %>% 
  dplyr::select(icgc_donor_id, sample_id, sample_id2)

TCGA_names_conversion_table = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/sample_ids_tables/TCGA_sampleid_conversion.tsv")
full_TCGA_names = data.frame(icgc_donor_id = as.character(),
                             sample_id2 = as.character())
for(TCGA_sample in TCGA_names_conversion_table$sample_id){
  
  # find full sample name in ../../data/
  full_sample_name = list.files('/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data', pattern = TCGA_sample)
  
  if(length(full_sample_name) >= 1){
    new_entry = data.frame(icgc_donor_id = as.character(filter(TCGA_names_conversion_table, sample_id == TCGA_sample) %>% pull(icgc_donor_id)),
                           sample_id2 = as.character(full_sample_name)) %>% 
      mutate(sample_id2 = gsub("muts_pass_", "", sample_id2),
             sample_id2 = gsub(".csv", "", sample_id2))
    # append to table
    full_TCGA_names = bind_rows(full_TCGA_names, new_entry)
  }
}
TCGA_names_conversion_table = merge(TCGA_names_conversion_table, full_TCGA_names) %>% 
  dplyr::select(icgc_donor_id, sample_id, sample_id2)

names_conversion_table = bind_rows(PCAWG_names_conversion_table, TCGA_names_conversion_table) %>% 
  mutate(sample_id = ifelse(!is.na(icgc_donor_id),
                            icgc_donor_id,
                            sample_id)) %>% 
  dplyr::select(sample_id, sample_id2)


## lymphoid samples' SBS84-like exposures obtained with sigprofiler
SBS84_exposures_lymph = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/AID_SHM/3_extract_SBS84__sort_AID_vs_nonAID_samples/lymphoid_samples_SBS84_exposure_K3.tsv") %>% 
  mutate(sample_id = sub('.*_', '', Sample)) %>% 
  select(sample_id, contains("SBS84")) %>% 
  rename("SBS84" = contains("SBS84"))
# non-lymphoid
SBS84_exposures_nonlymph = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/AID_SHM/5_non_lymphoid_controls/3_extract_SBS84__sort_AID_vs_nonAID_samples/nonlymphoid_samples_SBS84_exposure_K10.tsv") %>% 
  mutate(sample_id = sub('.*_', '', Sample)) %>% 
  select(sample_id, contains("SBS84")) %>% 
  rename("SBS84" = contains("SBS84"))
# bind them
SBS84_exposures = bind_rows(SBS84_exposures_lymph, SBS84_exposures_nonlymph)
  
## N A3A-putative TpC2DpH mut pairs
TpC2DpH_mutpairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv")


metadata = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv") %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x)) %>% 
  mutate(dataset = ifelse(str_detect(sample_id, "MSM0"),
                          "Kucab et al. 2019",
                          ifelse(str_detect(sample_id, "MSK0"),
                                 "Zou et al. 2021",
                                 source)),
         info1 = ifelse(str_detect(sample_id, "MSK0"),
                        paste0(info1, "ko"),
                        info1),
         info2 = gsub("^[a-z]_", "", info2), 
         info2 = ifelse(str_detect(sample_id, "MSM0") | str_detect(sample_id, "MSK0"),
                        gsub("DNA damage response inhibitors", "DNA damage resp. inh.", info2),
                        info2),
         info2 = ifelse(str_detect(sample_id, "MSK0"),
                        paste0(info2, " (pathway)"),
                        ifelse(str_detect(sample_id, "MSM0"),
                               paste0(info2, " (treatment)"),
                               info2))) %>% 
  rename("sample_id2" = "sample_id") %>% 
  merge(names_conversion_table, all = T) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            ifelse(!is.na(sample_id2),
                                   sample_id2,
                                   sample_id_2),
                            sample_id)) 

metadata_SBS84_exposures = metadata %>% 
  merge(SBS84_exposures, all = T) %>% 
  select(-c(sample_id, sample_id_2)) %>% 
  rename("sample_id" = "sample_id2") %>% 
  merge(SBS84_exposures, all = T) %>% 
  merge(names_conversion_table, all = T) %>% 
  select(-sample_id2) %>% 
  distinct
  
metadata_TpC2DpH_mutpairs = metadata %>% 
  merge(TpC2DpH_mutpairs, all = T) %>% 
  mutate(sample_id = ifelse(sample_id!=sample_id2 & !is.na(sample_id2),
                            sample_id2,
                            sample_id)) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarise_all(funs(toString(na.omit(.)))) %>% 
  mutate_all(na_if, "") %>% 
  distinct

metadata_SBS84_TpC2DpH = merge(metadata_SBS84_exposures,
                               metadata_TpC2DpH_mutpairs, all = T) %>% 
  filter(!(is.na(SBS84) & is.na(A3A_1kbpairs_TpC2DpH))) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarise_all(funs(toString(na.omit(.)))) %>% 
  mutate_all(na_if, "")

results_regressions = lapply(c("../../1_parser_and_regressions/res/results_real_samples.tsv",
                               "../../1_parser_and_regressions/res/simulated_positive_controls.tsv"),
                             read_tsv) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>% 
  rename_with(~str_replace(., 'high$', '')) %>% 
  rename_with(~str_replace(., 'hairpin_TpCpH$', '')) %>% 
  rename_with(~str_replace(., 'AID_target$', '')) %>%
  rename("sample_id2" = "sample_id") %>% 
  left_join(metadata_SBS84_TpC2DpH) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            sample_id2,
                            sample_id)) %>% 
  select(-c(info1,info2)) %>% 
  drop_na(starts_with("estimate_")) %>% 
  select(sample_id, contains("estimate_"), contains("conf"))



#### permutations

coefficient_table = results_regressions %>%
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat) %>% 
  group_by(sample_id, mark)


## Parameters and initializing of some objects
totalNumIters = 1000
coefficient_Resamp = list()
set.seed(1)

## Generate matrices resampling betas from their CI95% distributions (instead of UPmultinomial)
resample_from_CI = function(coefficient_table){
  summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                           mean = estimate,
                                                           sd = 1, 
                                                           lower = conf.low,
                                                           upper = conf.high))}
for (nIter in 1:totalNumIters) {
  print(paste0("nIter ", nIter))
  # for each sample (row) resample coefficients from CI95% distrs.
  coefficient_Resamp[[nIter]] = resample_from_CI(coefficient_table) %>% 
    mutate("nIter" = nIter)
  gc()
}

# bind all permutated tables as an input for the NMF
coefficient_Resamp = bind_rows(coefficient_Resamp)

write_tsv(coefficient_Resamp, paste0("permuted_coefficients_", totalNumIters, "iters.tsv"))




#####################################################################
#### NMF
#####################################################################


#### evaluate best combination of n variables and k signatures

coefficient_Resamp = coefficient_Resamp %>% 
  ungroup %>% 
  pivot_wider(names_from = mark, values_from = resampled_estimate) %>% 
  unite("sample_id", sample_id, nIter, sep = "__")

# a) keep positive coefficients, and convert negative to zero 
coefficient_Resamp_posmatrix = coefficient_Resamp %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))
# b) convert positive coefficients to zero, and convert negative to positive
coefficient_Resamp_negmatrix = coefficient_Resamp %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))
# merge converted coefficients into the NMF input
coefficient_Resamp = merge(coefficient_Resamp_posmatrix,
                           coefficient_Resamp_negmatrix,
                           by = "sample_id",
                           suffixes = c("_poscoeff", "_negcoeff")) %>%
  mutate(sample_id = gsub("__muts", "_muts", sample_id)) %>% 
  separate(sample_id, into = c("id", "nIter"), sep = "__", remove = F) %>% 
  mutate(sample_id = gsub("__", "_nIter", sample_id)) %>% 
  column_to_rownames("sample_id") %>% 
  select(-id) %>% 
  split(., f = .$nIter) %>% 
  map(.f = list(. %>% data.matrix))



## Run NMF for each matrix generated in the previous step

# Creates a set of copies of R running in parallel and communicating over sockets.
## Parameters and initializing of some objects
set.seed(1)
maxK = length(unique(coefficient_table$mark)) # max number of signatures to consider, it will go from 2 to maxK -- shouldn't be larger than nº of features
nCPUs = 8
cl = makeCluster(nCPUs)
clusterExport(cl = cl, list("coefficient_Resamp")) #, envir=tand)  # This can be slow for large lists
clusterEvalQ(cl, library(NMF))
# prepare the resampled NMF input matrices, generated just once and then re-used for any possible # factors/clusters
nmfHmatAllByFact = list()
nmfWmatAllByFact = list()

for (nFact in 2:maxK) {
  cat(sprintf("Running NMF: nFact %d (all iters)\n", nFact))
  nmfOutputByIter = parLapply(cl,
                              X = coefficient_Resamp,
                              fun = function(x, nFact) {
                                # store nIter
                                nIter = x[1,1]
                                x = x[,-1]
                                ## cannot handle all-zero columns (such as MSH6-neg, since all coeff. are +...)
                                # store their colnames
                                all_zero_columns = data.frame(x) %>% select_if(colSums(.) <= 0) %>% colnames
                                # remove them from matrix before NMF
                                x_nozerocols = x[, colSums(x) > 0]
                                ### NMF
                                nmf_run = nmf(x_nozerocols, rank = nFact, maxIter = 10000, seed = nIter)
                                # re-add the all-zero columns as zeros
                                all_zero_columns_matrix = data.frame(matrix(rep(0, len = length(all_zero_columns)), 
                                                                            ncol = length(all_zero_columns), 
                                                                            nrow = nrow(nmf_run@fit@H))) %>% 
                                  `colnames<-`(all_zero_columns) %>% as.matrix
                                full_matrix = cbind(nmf_run@fit@H, all_zero_columns_matrix) %>% 
                                  data.frame %>% 
                                  # reorder column names to be the same order as before
                                  select(all_of(colnames(coefficient_Resamp[[1]])[-1])) %>%
                                  as.matrix
                                nmf_run@fit@H = full_matrix
                                nmf_run},
                              nFact)
  idString = sprintf("nFact=%03d", nFact)
  nmfHmatAllByFact[[idString]] = matrix(nrow = 0, ncol = ncol(nmfOutputByIter[[1]])); # columns = original features...
  nmfWmatAllByFact[[idString]] = matrix(nrow = 0, ncol = nrow(nmfOutputByIter[[1]])); # columns = samples...
  for (nIter in 1:totalNumIters) {
    nmfOutput = nmfOutputByIter[[nIter]]
    rownames(nmfOutput@fit@H) = sprintf("run%03d_fact%02d", nIter, 1:nFact)
    colnames(nmfOutput@fit@W) = sprintf("run%03d_fact%02d", nIter, 1:nFact)
    nmfHmatAllByFact[[idString]] = nmfHmatAllByFact[[idString]] %>% rbind(nmfOutput@fit@H) # same thing that you get using the coef(NMFfit) command in R   -> loadings
    nmfWmatAllByFact[[idString]] = nmfWmatAllByFact[[idString]] %>% rbind(t(nmfOutput@fit@W)) # same thing that you get using the basis(NMFfit) command in R  -> transformed data
  }
}
stopCluster(cl)
rm(idString, nmfOutput, nFact)


## Subtract results of pos - results of neg, and get the absolute
nmfHmatAllByFact_combined_pos_negcoeff = nmfHmatAllByFact %>% 
  map(., ~as_tibble(.x, rownames = NA)) %>%
  map(., ~rownames_to_column(.x, "id")) %>%
  map(., ~pivot_longer(.x,
                       cols = contains("coeff"),
                       names_to = "dna_repair_mark",
                       values_to = "Weight")) %>%
  map(., ~mutate(.x, dna_repair_mark = gsub("_...coeff", "", dna_repair_mark))) %>% 
  map(., ~group_by(.x, dna_repair_mark, id)) %>%
  map(., ~summarise(.x, Weight = abs(.Primitive("-")(Weight[1], Weight[2])))) %>%
  map(., ~ungroup(.x)) %>%
  ## scale weights per signature so they add up to 1
  map(., ~group_by(.x, id)) %>% 
  map(., ~mutate(.x, sumWeight = sum(Weight))) %>% 
  map(., ~group_by(.x, id, dna_repair_mark)) %>% 
  map(., ~summarise(.x, Weight = Weight/sumWeight)) %>% 
  map(., ~ungroup(.x)) %>% 
  ## go to original format
  map(., ~pivot_wider(.x, names_from = dna_repair_mark, values_from = Weight)) %>%
  # arrange iterations (rows)
  map(., ~arrange(.x, id)) %>%
  map(., ~column_to_rownames(.x, 'id')) %>%
  # arrange dna repair mark names (columns)
  map(.f = list(. %>% select(coefficient_Resamp[[1]] %>% 
                               data.frame %>% 
                               select(contains("poscoeff")) %>% 
                               names() %>% 
                               gsub("_poscoeff", "", .)))) %>%
  map(., ~as.matrix(.x))



## Cluster the NMF factors, and find combination that yields best silhouette (and possibly deviation from random data)
## Also run in parallel on as many CPUs as specified
## Returns a data.table with nFact, k, and the quality measures per nFact*k combination.
#' @param maxNumClusters Global override - holds true for any number of nFact... note that some nFact will have less clusters than that 
#'                       For a given nFact, the max # clusters is 5*nFact.  Meaning: max 10 clusters for 2 factors, max 15 for 3 factors, etc.
#'                       (probably even this many does not make sense to extract...)

#initialize cosine distance function which is used for clustering
cosineDist <- function(x) {
  mat = as.dist(t(1 - cosine(t(x))))
  return(mat)
}

#initialize function for clustering
optimizeClustersPar = function(nmfHmatAllByFact = NULL, maxNumClusters = maxK, threads = nCPUs) {
  # the highest nFact (number of factors) we had checked in this dataset... the lowest is always assumed to be ==2 but this could be inspected as well
  maxNumFact = names(nmfHmatAllByFact) %>% str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer %>% max
  
  cl = makeCluster(threads, useXDR = FALSE);
  clusterEvalQ(cl, library(cluster)) # to be able to run "pam"  (the library 'cluster' refers to data clustering, not compute-clusters!)  
  sigDiscoScores = data.table()
  numNmfIters = -1; # this will be found from the nmfHmatAllByFact
  
  for (aNmfSettings in names(nmfHmatAllByFact)) {
    
    # aNmfSettings is a string like "nFact=002", which encodes both the "nFact" parameter
    nFact = aNmfSettings %>% str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer
    
    # how many iterations of NMF were run for this maxNumFact.
    if (numNmfIters == -1) {
      numNmfIters = nmfHmatAllByFact[[aNmfSettings]] %>% rownames %>% str_match("^run(\\d+)_") %>% .[, 2] %>% as.integer %>% max
    } else {
      numNmfIters2 = nmfHmatAllByFact[[aNmfSettings]] %>% rownames %>% str_match("^run(\\d+)_") %>% .[, 2] %>% as.integer %>% max
      if (numNmfIters != numNmfIters2) {
        cat(sprintf("For aNmfSettings=%s, the numNmfIters %d is different than previously, when it was %d. Will continue running normally.\n", aNmfSettings, numNmfIters2, numNmfIters))
      }
    }
    
    nmfReBootstrappedDists = cosineDist(nmfHmatAllByFact[[aNmfSettings]])
    
    # parallel code: in each thread, runs the "pam" clustering for a different # clusters, starting from the same matrix of cosine similarities
    # (note: Sanger (Alexandrov et al.) always run this with #clusters = #NMF factors, which may not be ideal)
    clusterExport(cl = cl, list("nmfReBootstrappedDists"), envir = environment()) # , envir=tand)
    pamOutputByIter = parLapplyLB(# load balancing - use parLapplyLB
      cl = cl,
      X = as.list(rev(2:min(maxNumClusters, nFact * 5))), # pam() with different values of k is run in parallel 
      fun = function(x) {
        set.seed(42 + x);
        pam(nmfReBootstrappedDists, k = x, diss = T)
      })
    
    for (clusters in pamOutputByIter) {
      k = length(clusters$medoids);
      sigDiscoScores = sigDiscoScores %>% rbind(data.table(nFact = nFact, k = k,
                                                           avgClScore = mean(clusters$silinfo$clus.avg.widths), minClScore = min(clusters$silinfo$clus.avg.widths),
                                                           secondWorstClScore = clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                                                           avgPtScore = mean(silhouette(clusters)[, "sil_width"]), medianPtScore = median(silhouette(clusters)[, "sil_width"])
      ))
      cat(
        sprintf("%02d\t%02d\t%.4f\t%.4f\t%.4f\t%.4f\n",
                nFact, k,
                mean(clusters$silinfo$clus.avg.widths), min(clusters$silinfo$clus.avg.widths),
                clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                mean(silhouette(clusters)[, "sil_width"]), median(silhouette(clusters)[, "sil_width"])
        )
      )
    }
  }
  
  stopCluster(cl)
  #rm(maxNumClusters, threads, cl, pamOutputByIter, k )
  return(sigDiscoScores)
}

# run clustering
sigClusterScores = optimizeClustersPar(nmfHmatAllByFact_combined_pos_negcoeff, maxNumClusters = maxK)

## Visualize the clustering quality scores in heatmaps (facetted by the 'smooth' parameter)
heatmap_clustering = sigClusterScores[k <= maxK] %>% 
  melt(id.vars = c("nFact", "k"),
       measure.vars = "avgClScore") %>%
  rename("stability" = "value") %>% 
  ggplot(aes(nFact, k)) +
  geom_tile(aes(fill = stability)) +
  geom_text(aes(label = round(stability, 2))) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.4)
ggsave("NMF_heatmap_clustering.jpg",
       plot = heatmap_clustering,
       device = "jpg",
       width = 12.5,
       height = 6,
       dpi = 600)



#################################################################################
###### now run the final NMF knowing the optimal k and n features


##### split the original coefficients between pos and neg

# a) keep positive coefficients, and convert negative to zero 
coefficients_posmatrix = coefficient_table %>% 
  select(sample_id, mark, contains("estimate")) %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))
# b) convert positive coefficients to zero, and convert negative to positive
coefficients_negmatrix = coefficient_table %>% 
  select(sample_id, mark, contains("estimate")) %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))

# regenerate the coeff matrix and transpose, for RcppML::nmf()
coefficient_matrix_RcppML = bind_rows(mutate(coefficients_posmatrix, submatrix = "poscoeff"),
                                      mutate(coefficients_negmatrix, submatrix = "negcoeff")) %>% 
  unite("dna_repair_mark",mark, submatrix, sep = "_") %>% 
  mutate(dna_repair_mark = gsub("estimate_", "", dna_repair_mark)) %>% 
  # transpose
  pivot_wider(names_from = sample_id, values_from = estimate) %>% 
  column_to_rownames("dna_repair_mark") %>% 
  as.matrix


##### prepare condition_pathway_pairs for "good-model-score"
repair_mark_pathways = coefficient_table %>% 
  select(sample_id,  mark, estimate) %>% 
  rename("dna_repair_mark" = "mark",
         "sample_id" = "sample_id") %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(dna_repair_mark, "OGG1_"),
                                                   "BER",
                                                   ifelse(str_detect(dna_repair_mark, "UV_"),
                                                          "NER",
                                                          ifelse(str_detect(dna_repair_mark, "MSH6_|SETD2_"),
                                                                 "MMR",
                                                                 ifelse(str_detect(dna_repair_mark, "XRCC4"),
                                                                        "DSBR",
                                                                        "other"))))) # ,#ADD TP53
condition_pathway_pairs = metadata_SBS84_TpC2DpH %>% 
  select(info1, info2) %>% 
  distinct %>% 
  filter(! str_detect(info2, "[C,c]ontrol")) %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(info1, "[G,g]amma"),
                                                   "BER,DSBR",
                                                   ifelse(str_detect(info2, "BER |[A,a]lkylating|[N,n]itrosamine"),
                                                          "BER",
                                                          ifelse(str_detect(info2, "NER |[A,a]romatic|[H,h]eterocyclic|PAH|Radiation") |
                                                                   str_detect(info1, "platin"),
                                                                 "NER",
                                                                 ifelse(str_detect(info2, "MMR "),
                                                                        "MMR",
                                                                        ifelse(str_detect(info2, "DSB|DNA damage resp. inh.|[H,h]elicas|NHEJ|MMEJ|HR ") |
                                                                                 str_detect(info1, "PARP|Etoposide|Bleomycin|Camptothecin|Olaparib|Temozolomide|Melphalan|Cyclophosphamide|Mechlorethamine"),
                                                                               "DSBR",
                                                                               "other")))))) %>% # ADD TP53
  separate_rows(putative_repair_pathway_involved, sep = ",") %>% 
  filter(!is.na(putative_repair_pathway_involved)) %>% 
  left_join(metadata_SBS84_TpC2DpH) %>% 
  filter(sample_id %in% coefficient_table$sample_id) %>% 
  arrange(putative_repair_pathway_involved) %>% 
  left_join(repair_mark_pathways) %>% 
  select(sample_id, dna_repair_mark) %>% 
  rename("DNA repair\nactivity" = "dna_repair_mark") %>% 
  group_by(sample_id) %>%
  mutate(sample_id_possible_marks = n()) %>%
  ungroup()
condition_pathway_pairs = condition_pathway_pairs %>% 
  rowwise %>% 
  mutate(max_sample_mark_score = 1 / length(unique(condition_pathway_pairs$sample_id)) / sample_id_possible_marks) %>%
  select(-c(sample_id_possible_marks))


####
### Run NMF for several k's

# for plots
jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))

for(optimal_k in seq(5, maxK)){
  
  # final NMF (here using RcppML::nmf instead of NMF::nmf as above)
  nmf_res = RcppML::nmf(coefficient_matrix_RcppML, 
                        k = optimal_k,
                        maxit = 10000, 
                        seed = 1)
  
  # Signature exposures in samples
  exposures = nmf_res$h %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column("Signature") %>%
    pivot_longer(cols = !contains("Signature"), names_to = "sample_id", values_to = "Exposure") %>% 
    # add metadata info (e.g. treatments, MSI, HR, smoking...)
    merge(metadata_SBS84_TpC2DpH, all = T) %>% 
    filter(!is.na(Exposure)) %>% 
    separate(sample_id, into = c("sample_id", "mut x-fold"), sep = "__", fill = "right") %>% 
    # simulated pos controls dont have SBS84 nor A3A mut pairs, so add them an average fake one so that their size is not tiny in the plots (as point size is dependent on SBS84 raw exposure)
    mutate(SBS84 = ifelse(is.na(SBS84) & !is.na(`mut x-fold`),
                          mean(metadata$SBS84, na.rm = T),
                          SBS84)) %>%
    mutate(A3A_1kbpairs_TpC2DpH = ifelse(is.na(A3A_1kbpairs_TpC2DpH) & !is.na(`mut x-fold`),
                                        mean(metadata$A3A_1kbpairs_TpC2DpH, na.rm = T),
                                        A3A_1kbpairs_TpC2DpH)) %>%
    mutate(sample_id = gsub("-high_activity", "", sample_id),
           dataset = ifelse(is.na(dataset),
                            `mut x-fold`,
                            dataset),
           Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h))))) %>%
    # highlight top hits for each signature
    group_by(Signature) %>% 
    mutate(is.hit = ifelse(Exposure==max(Exposure), "Top hit", NA))
  
  # DNA repair mark weights in signatures
  weights = nmf_res$w %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column("dna_repair_mark") %>%
    pivot_longer(cols = contains("nmf"), names_to = "signature", values_to = "weight") %>% 
    extract(dna_repair_mark, into = c("dna_repair_mark", "submatrix"), "(.*)_([^_]+$)") %>% 
    mutate(dna_repair_mark = factor(dna_repair_mark, levels = unique(gsub("_...coeff", "", rownames(nmf_res$w)))),
           signature = factor(signature, levels = unique(exposures$Signature))) %>% 
    arrange(dna_repair_mark, signature) %>% 
    group_by(dna_repair_mark, signature) %>% 
    ## subtract results of pos - results of neg, and get the absolute
    summarise(Weight = abs(.Primitive("-")(weight[1], weight[2]))) %>% 
    ungroup %>% 
    rename("DNA repair\nactivity" = "dna_repair_mark", 
           "Signature" = "signature") %>% 
    relocate(Signature)
  
  #####
  ## calculate "good-model-score": 1 would mean that ALL samples have 100% exposure from the signature(s) that is contributed 100% by a specific DNArep mark whose mechanism/pathway overlaps completely with the sample´s condition (gene-/-, treatment...), and 0 the opposite
  # actually only considering "ground-truth" sample_condition--DNArep_pathway pairs ('condition_pathway_pairs')
  good_model_score_table = weights %>% 
    merge(condition_pathway_pairs, all = T) %>% 
    merge(exposures, all = T) %>% 
    rename("sensical sample-mark pair" = "DNA repair\nactivity",
           "mark weight" = "Weight",
           "signature exposure" = "Exposure") %>% 
    select(sample_id, "sensical sample-mark pair", Signature, "mark weight", "signature exposure", max_sample_mark_score) %>% 
    mutate(sample_mark_score = `mark weight` * `signature exposure` * max_sample_mark_score) %>% 
    arrange(sample_id, `sensical sample-mark pair`, Signature)
  
  #write_tsv(good_model_score_table, paste0("K", optimal_k, "_table.tsv"))
  
  good_model_score = drop_na(good_model_score_table) %>% pull(sample_mark_score) %>% sum
  #####
  
  
  ### plotting
  
  # for each nmf signature, plot only the top n samples with highest Exposure
  top_n_samples = 5
  filtered_exposures = exposures %>% 
    group_by(Signature) %>% 
    arrange(desc(Exposure)) %>% 
    slice_head(n = top_n_samples) %>% 
    ungroup
  
  ## exposures (only top_n_samples)
  pos = position_jitter(w = 0.25, h = 0, seed = 1)
  exposures_plot = ggplot(filtered_exposures %>% 
                            mutate(`SBS84-like exposure` = as.numeric(SBS84),
                                   `5'-T(C>D)H-3' pairs` = as.numeric(A3A_1kbpairs_TpC2DpH)), 
                          aes(x = Signature,
                              y = Exposure*100)) +
    scale_y_continuous(labels = function(x) sub("0+$", "", x)) +
    geom_point(aes(fill = `SBS84-like exposure`,
                   size = `5'-T(C>D)H-3' pairs`),
                   shape = 21, #dataset, ##scale_shape_manual(values = c(21, 23, 25, seq(1,length(unique(filtered_exposures$dataset))-3,1))) +
                   position = pos) +
    scale_fill_gradient(low = "white", high = "black") +
    guides(fill = guide_legend(override.aes = list(size=4, shape=21)),
           shape = guide_legend(override.aes = list(size=4))) +
    ggrepel::geom_text_repel(aes(label = sample_id),
                             size = 3,
                             force = 5,
                             position = pos,
                             max.overlaps = 1000000,
                             min.segment.length = 1) +
    facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
    theme_classic() +
    xlab("") +
    #ggtitle(paste0("Model Score = ", round(good_model_score*100, 2), "%")) +
    ylab(paste0("% Exposure (top-", top_n_samples, " samples)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(6, "mm"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 8))
  
  ## weights
  weights_plot = ggplot(weights %>%
                          mutate(`DNA repair\nactivity` = gsub("_2strands", "", `DNA repair\nactivity`),
                                 Signature = factor(gsub("nmf", "", Signature), levels = seq(1:length(levels(weights$Signature))))) %>% 
                          # scale weights per signature so they add up to 1
                          group_by(Signature) %>% 
                          mutate(sumWeight = sum(Weight)) %>% 
                          group_by(Signature, `DNA repair\nactivity`) %>% 
                          summarise(Weight = Weight/sumWeight) %>% 
                          ungroup, 
                        aes(x = Signature,
                            y = Weight)) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) sub("0+$", "", x)) +
    geom_col(aes(fill = `DNA repair\nactivity`)) +
    scale_fill_manual(values = jet.colors(length(levels(weights$`DNA repair\nactivity`)))) +
    guides(fill = guide_legend(override.aes = list(size=6))) +
    facet_grid(cols = vars(Signature), scales = "free") +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(2, "mm"),
          legend.text = element_text(size = 8))
  
  combined_plots = cowplot::plot_grid(NULL,
                                      cowplot::plot_grid(exposures_plot, NULL, nrow = 1, rel_widths = c(1, 0.04*(optimal_k/11))),
                                      NULL,
                                      weights_plot,
                                      nrow = 4,
                                      rel_heights = c(0.02, 0.75,-0.05,1))
  ggsave(paste0("plots/NMF_exposures_weights_plot_k", optimal_k, ".jpg"),
         plot = combined_plots,
         device = "jpg",
         width = 21.3,
         height = 12,
         dpi = 700,
         bg = "white")
}
