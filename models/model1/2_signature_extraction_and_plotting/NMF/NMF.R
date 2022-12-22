library(tidyverse)
library(data.table)
library(sampling)
library(parallel)
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


dir.create("plots")


##### parse input

metadata = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv") %>% 
  rename("Sample" = "sample_id") %>% 
  mutate(info2 = gsub("DNA damage response inhibitors",
                               "DNA damage resp. inh.",
                               info2))

# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../1_parser_and_regressions/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_"))

# a) keep positive coefficients, and convert negative to zero 
results_regressions_posmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))

# b) convert positive coefficients to zero, and convert negative to positive
results_regressions_negmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))


#### evaluate best combination of n variables and k signatures

# merge converted coefficients into the NMF input
coefficient_matrix = merge(results_regressions_posmatrix,
                           results_regressions_negmatrix,
                           by = "sample_id",
                           suffixes = c("_poscoeff", "_negcoeff")) %>%
  column_to_rownames("sample_id") %>% 
  rename_all(~str_replace_all(., 'estimate_', '')) %>% 
  ## WARNING: I'm MULTIPLYING all coefficients by 100, so there are no decimals (UPmultinomial only works with whole numbers, so it converts e.g. 0.2 to 0...)
  mutate(across(where(is.numeric), 
                ~.x * 100)) %>% 
  as.matrix

## Parameters and initializing of some objects
totalNumIters = 100
maxK = length(colnames(coefficient_matrix)) / 2  # max number of signatures to consider, it will go from 2 to maxK -- shouldn't be larger than nº of features, which is length(colnames(coefficient_matrix)) / 2 because matrix is duplicated into 2 submatrices (positive and negative coefficients)
nCPUs = 8
set.seed(1)

# prepare the resampled NMF input matrices, generated just once and then re-used for any possible # factors/clusters
nmfHmatAllByFact = list()
nmfWmatAllByFact = list()
coefficient_matrix_Resamp = list()


## Generate matrices with small random perturbations of the original matrix 
for (nIter in 1:totalNumIters) {
  cat(sprintf("Generating bootstrap matrix: nIter %d\n", nIter))
  coefficient_matrix_TempIter = matrix(data = NA,
                                       nrow = nrow(coefficient_matrix), ncol = ncol(coefficient_matrix),
                                       dimnames = list(rownames(coefficient_matrix), colnames(coefficient_matrix))
  )
  for (aRow in 1:nrow(coefficient_matrix)) {
    # for each sample (row in input matrix) change coefficients a bit, while keeping the total sum the same
    coefficient_matrix_TempIter[aRow,] = UPmultinomial(coefficient_matrix[aRow,]) / 100 # WARNING: dividing back by 100 to compensate the *100 in the coefficient_matrix making step
  }
  attr(coefficient_matrix_TempIter, "nIter") = nIter
  coefficient_matrix_Resamp[[nIter]] = coefficient_matrix_TempIter;
}
rm(coefficient_matrix_TempIter, nIter)


## Run NMF for each matrix generated in the previous step

# Creates a set of copies of R running in parallel and communicating over sockets.
cl = makeCluster(nCPUs)
clusterExport(cl = cl, list("coefficient_matrix_Resamp")) #, envir=tand)  # This can be slow for large lists
clusterEvalQ(cl, library(NMF))

for (nFact in 2:maxK) {
  cat(sprintf("Running NMF: nFact %d (all iters)\n", nFact))
  nmfOutputByIter = parLapply(cl,
                              X = coefficient_matrix_Resamp,
                              fun = function(x, nFact) {
                                   ## cannot handle all-zero columns (such as MSH6-neg, since all coeff. are +...)
                                   # store their colnames
                                   all_zero_columns = colnames(x[, colSums(x) <= 0])
                                   # remove them from matrix before NMF
                                   x_nozerocols = x[, colSums(x) > 0]
                                   # NMF
                                   nmf_run = nmf(x_nozerocols, rank = nFact, maxIter = 10000, seed = attr(x, "nIter"))
                                   # re-add the all-zero columns as zeros
                                   all_zero_columns_matrix = data.frame(matrix(rep(0, len = length(all_zero_columns)), 
                                                                               ncol = length(all_zero_columns), 
                                                                               nrow = nrow(nmf_run@fit@H))) %>% 
                                     `colnames<-`(all_zero_columns) %>% as.matrix
                                   full_matrix = cbind(nmf_run@fit@H, all_zero_columns_matrix) %>% data.frame %>% 
                                     # reorder column names to be the same order as before
                                     select(all_of(colnames(coefficient_matrix_Resamp[[1]]))) %>% as.matrix
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


## NEW: subtract results of pos - results of neg, and get the absolute
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
  map(.f = list(. %>% select(results_regressions_posmatrix %>% 
                               select(contains("estimate")) %>% 
                               names() %>% 
                               gsub("estimate_", "", .)))) %>%
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
ggsave("plots/NMF_heatmap_clustering.jpg",
       plot = heatmap_clustering,
       device = "jpg",
       width = 12.5,
       height = 6,
       dpi = 600)



###### now run the final NMF using the optimal k and n features

# regenerate the coeff matrix (without mult. by 100) and transpose, for RcppML::nmf()
coefficient_matrix_RcppML = bind_rows(mutate(results_regressions_posmatrix, submatrix = "poscoeff"),
                                      mutate(results_regressions_negmatrix, submatrix = "negcoeff")) %>% 
  pivot_longer(cols = contains("estimate"),
               names_to = "dna_repair_mark",
               values_to = "estimate") %>% 
  unite("dna_repair_mark", dna_repair_mark, submatrix, sep = "_") %>% 
  mutate(dna_repair_mark = gsub("estimate_", "", dna_repair_mark)) %>% 
  # transpose
  pivot_wider(names_from = sample_id, values_from = estimate) %>% 
  column_to_rownames("dna_repair_mark") %>% 
  as.matrix

for(optimal_k in seq(2, maxK)){
  
  # final NMF (here using RcppML::nmf instead of NMF::nmf as above)
  nmf_res = RcppML::nmf(coefficient_matrix_RcppML, 
                        k = optimal_k,
                        maxit = 10000, 
                        seed = 1)
  
  # Signature exposures in samples
  exposures = nmf_res$h %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column("Signature") %>%
    pivot_longer(cols = !contains("Signature"), names_to = "Sample", values_to = "Exposure") %>% 
    # add metadata info (e.g. treatments, MSI, HR, smoking...)
    left_join(metadata) %>% 
    mutate(Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h))))
           #MSI_parsed = ifelse(MSI_status %in% c("HYPER", "MSI", "ERCC2mut"), "MSI", NA),
           #hr_parsed = ifelse(hr_status %in% c("HR_deficient"), "HRdef", NA),
           #smoking_history_parsed = ifelse(smoking_history %in% c("Current", "Former"), "Smoker", NA),
           #treatment_platinum_parsed = ifelse(treatment_platinum == "TRUE", "Platinum", NA),
           #treatment_5FU_parsed = ifelse(treatment_5FU == "TRUE", "5FU", NA)
           ) %>%
    #unite(col = "Metadata", MSI_parsed, hr_parsed, smoking_history_parsed,treatment_platinum_parsed, treatment_5FU_parsed, na.rm = T, sep = " & ") %>% 
    #mutate(Metadata = gsub("^$", "NOTA/NA", Metadata)) %>% 
    #rename("Database" = "source") %>% 
    # highlight top hits for each signature
    group_by(Signature) %>% 
    mutate(is.hit = ifelse(Exposure==max(Exposure), "hit", NA)
           #has.Metadata = ifelse(Metadata != "NOTA/NA", "yes", "no")
           ) %>% 
    rename("Treatment\ntype" = "info2")
  # # "NOTA/NA" to the end
  # exposures$Metadata = factor(exposures$Metadata, levels = c(unique(exposures$Metadata)[unique(exposures$Metadata) != "NOTA/NA"],
  #                                                            "NOTA/NA"))
  # write_tsv(exposures,
  #           "NMF_exposures.tsv")
  
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
  # write_tsv(weights,
  #           "NMF_weights.tsv")
  
  
  ### plotting
  
  #jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))
  
  exposure_limit = 0.025
  
  # plot only samples with Exposures > exposure_limit%
  filtered_exposures = filter(exposures, Exposure > exposure_limit)
  
  # detect signatures with no sample having at least the exposure_limit exposure to it, to show it in the plot
  signatures_lower_exposure_than_limit = exposures %>% 
    group_by(Signature) %>% 
    summarise(max_exposure = max(Exposure)) %>% 
    filter(max_exposure < exposure_limit) %>% pull(Signature) %>% as.character
  
  if(length(signatures_lower_exposure_than_limit) > 0){
    filtered_exposures = filtered_exposures %>% 
      bind_rows(data.frame(signatures_lower_exposure_than_limit) %>% 
                  `colnames<-`("Signature"))
  }
  
  ## exposures (only samples with Exposures > exposure_limit%)
  exposures_plot = ggplot(filtered_exposures, 
                          aes(x = Signature,
                              # convert exposures to %
                              y = Exposure*100,
                              group = `Treatment\ntype`)) + # Database
    scale_y_continuous(labels = function(x) sub("0+$", "", x)) +
    # # all points, no colors
    # geom_point(aes(shape = Database),
    #            position = position_dodge(width = 1),
    #            alpha = 0.5,
    #            size = 4) +
    # # colored those with metadata info
    geom_point(aes(#shape = Database,
                   fill = `Treatment\ntype`
                   #alpha = has.Metadata
                   ),
               shape = 21,
               #position = position_dodge(width = 1),
               size = 4) +
    #scale_shape_manual(values = c(21,23,24,22,20,25)) + #"\u2716"
    scale_fill_manual(values = jet.colors(length(unique(exposures$`Treatment\ntype`)))) +
    #scale_fill_manual(values = c(jet.colors(length(unique(filter(exposures, Exposure > 0.001 & Metadata!="NOTA/NA")$Metadata))), "white")) + # Metadata==NOTA/NA is assigned white
    #scale_alpha_manual(values = c(0, 0.8), guide = 'none') +
    guides(shape = guide_legend(override.aes = list(size=6)),
           fill = guide_legend(override.aes = list(size=6, shape=21))) +
    # label sample names
    ggrepel::geom_text_repel(data = filtered_exposures %>%  #filter(is.hit == "hit") %>%
                                    mutate(Sample = gsub("MSM0.", "samp_", Sample)),
                             aes(label = Sample),
                             size = 4,
                             #nudge_y = 5,
                             force = 5,
                             max.overlaps = 1000000,
                             min.segment.length = 10000) +
    facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
    theme_classic() +
    xlab("") +
    ylab(paste0("% Exposure (>", exposure_limit*100, "%)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(4, "mm"),
          legend.text = element_text(size = 10))

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
          legend.text = element_text(size = 9))
  
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
