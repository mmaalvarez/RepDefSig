library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(cowplot)
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

dir.create("plots")

jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))


##### parse input

metadata = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info"))) %>% 
  rename("Sample" = "sample_id") %>% 
  mutate(dataset = ifelse(str_detect(Sample, "MSM0"),
                          "Kucab et al. 2019",
                          ifelse(str_detect(Sample, "MSK0"),
                                 "Zou et al. 2021",
                                 "ERROR: Unexpected sample name")),
         info1 = ifelse(str_detect(Sample, "MSK0"),
                        paste0(info1, "ko"),
                        ifelse(str_detect(Sample, "MSM0"),
                               info1,
                               "ERROR: Unexpected sample name")),
         info2 = gsub("^[a-z]_", "", info2), 
         info2 = gsub("DNA damage response inhibitors", "DNA damage resp. inh.", info2),
         info2 = ifelse(str_detect(Sample, "MSK0"),
                        paste0(info2, " (pathway)"),
                        ifelse(str_detect(Sample, "MSM0"),
                               paste0(info2, " (treatment)"),
                               "ERROR: Unexpected sample name")))

# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../1_parser_and_regressions/res/results.tsv") %>% #../1_parser_and_regressions/bin/output/results.tsv") %>% #../../model1/1_parser_and_regressions/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_"), contains("conf"), glm)



#########################################################################
### PCA to look for biases regarding which samples were regressed with which glm family
pca_res = results_regressions %>% 
  select(sample_id, contains("estimate_")) %>% 
  column_to_rownames("sample_id") %>% 
  rename_with(~str_replace(., 'estimate_', '')) %>% 
  rename_with(~str_replace(., '.L', '')) %>%
  PCA(graph = FALSE)

pca = results_regressions %>%
  merge(pca_res$ind$coord %>% 
          data.frame %>% 
          rownames_to_column("sample_id")) %>% 
  select(sample_id, contains("Dim"), glm) %>% 
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  rename("Sample" = "sample_id") %>% 
  merge(metadata) %>% 
  ggplot(aes(x = PC1,
             y = PC2)) +
  coord_fixed() +
  geom_point(aes(fill = info2,
                 shape = glm)) +
  stat_ellipse(geom = "polygon",
               aes(col = glm),
               alpha = 0,
               show.legend = T,
               level = 0.95) +
  scale_fill_manual(values = jet.colors(length(unique(metadata$info2)))) +
  scale_color_manual(values = jet.colors(length(unique(metadata$dataset)))) +
  scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(size=6, shape=21))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  ggrepel::geom_text_repel(aes(label = info1),
                           force = 5,
                           segment.size = 0.1,
                           min.segment.length = 0.001,
                           max.overlaps = 100000,
                           max.iter = 100000) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size=15),
        legend.title = element_blank())
ggsave("plots/pca.svg",
       plot = pca,
       device = "svg",
       width = 20,
       height = 15.6,
       dpi = 600,
       bg = "white")

scree = fviz_eig(pca_res,
                 ylim = c(0, round(max(pca_res$eig[,2]))),
                 geom = c("bar"),
                 addlabels = FALSE,
                 ncp = length(rownames(pca_res$eig)),
                 main = "",
                 ggtheme = theme_classic(base_size = 20),
                 xlab = "PC",
                 ylab = "% variance")
ggsave("plots/scree.svg",
       plot = scree,
       device = "svg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")

vars = pca_res$var$coord %>%
  data.frame() %>%
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  rownames_to_column("vars") %>% 
  ggplot() +
  geom_segment(aes(x = 0, xend = PC1,
                   y = 0, yend = PC2),
               arrow = arrow(length = unit(0.025,
                                           "npc"),
                             type = "open"),
               lwd = 0.5,
               linetype = "dashed") +
  ggrepel::geom_text_repel(aes(x = PC1,
                               y = PC2,
                               label = vars),
                           size = 6,
                           direction = "y",
                           vjust = 3,
                           force = 5,
                           segment.size = 0,
                           min.segment.length = 0,
                           max.iter = 100000) +
  theme_nothing()
ggsave("plots/vars.svg",
       plot = vars,
       device = "svg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")



#####################################################################
#### NMF
#####################################################################


###################################################################
#### evaluate best combination of n variables and k signatures

coefficient_table = results_regressions %>%
  select(sample_id, contains("estimate_"), contains("conf")) %>% 
  rename_all(~str_replace_all(., '.L', '')) %>% 
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat) %>% 
  group_by(sample_id, mark)

## Parameters and initializing of some objects
totalNumIters = 100
maxK = length(unique(coefficient_table$mark)) # max number of signatures to consider, it will go from 2 to maxK -- shouldn't be larger than nº of features
nCPUs = 8
set.seed(1)

# prepare the resampled NMF input matrices, generated just once and then re-used for any possible # factors/clusters
nmfHmatAllByFact = list()
nmfWmatAllByFact = list()
coefficient_matrix_Resamp = list()

## Generate matrices resampling betas from their CI95% distributions (instead of UPmultinomial)
resample_from_CI = function(coefficient_table){
  coefficient_table %>% 
    summarise(resampled_estimate = rtnorm(n = 1,
                                          mean = estimate,
                                          sd = 1, 
                                          lower = conf.low,
                                          upper = conf.high))
}

for (nIter in 1:totalNumIters) {
  cat(sprintf("Generating bootstrap matrix: nIter %d\n", nIter))
  
  # for each sample (row) resample coefficients from CI95% distrs.
  coefficient_matrix_TempIter = resample_from_CI(coefficient_table) %>% 
    pivot_wider(names_from = mark, values_from = resampled_estimate) %>% 
    ungroup
    
  # a) keep positive coefficients, and convert negative to zero 
  coefficient_matrix_TempIter_posmatrix = coefficient_matrix_TempIter %>% 
    mutate_if(is.numeric,
              ~if_else(.<=0, 0, .))
  # b) convert positive coefficients to zero, and convert negative to positive
  coefficient_matrix_TempIter_negmatrix = coefficient_matrix_TempIter %>% 
      mutate_if(is.numeric,
                ~if_else(.>=0, 0, abs(.)))
  # merge converted coefficients into the NMF input
  coefficient_matrix_TempIter = merge(coefficient_matrix_TempIter_posmatrix,
                                      coefficient_matrix_TempIter_negmatrix,
                                      by = "sample_id",
                                      suffixes = c("_poscoeff", "_negcoeff")) %>%
    column_to_rownames("sample_id") %>%
    data.matrix
  
  attr(coefficient_matrix_TempIter, "nIter") = nIter
  
  coefficient_matrix_Resamp[[nIter]] = coefficient_matrix_TempIter
  gc()
}
rm(coefficient_matrix_TempIter)


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
  map(.f = list(. %>% select(results_regressions_posmatrix %>% 
                             select(contains("estimate")) %>% 
                             names() %>% 
                             gsub("estimate_", "", .) %>% 
                             gsub(".L", "", .)))) %>%
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



#################################################################################
###### now run the final NMF knowing the optimal k and n features

## split the original results_regressions between pos and neg

# a) keep positive coefficients, and convert negative to zero 
results_regressions_posmatrix = results_regressions %>% 
  select(sample_id, contains("estimate_")) %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))
# b) convert positive coefficients to zero, and convert negative to positive
results_regressions_negmatrix = results_regressions %>% 
  select(sample_id, contains("estimate_")) %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))

# regenerate the coeff matrix and transpose, for RcppML::nmf()
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


## prepare condition_pathway_pairs for "good-model-score"
repair_mark_pathways = results_regressions %>% 
  select(contains("estimate")) %>% 
  rename_with(~str_replace(., 'estimate_', '')) %>% 
  rename_with(~str_replace(., '.L', '')) %>% 
  colnames %>% 
  data.frame %>% 
  `colnames<-`("dna_repair_mark") %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(dna_repair_mark, "OGG1_"),
                                                  "BER",
                                                  ifelse(str_detect(dna_repair_mark, "UV_"),
                                                         "NER",
                                                         ifelse(str_detect(dna_repair_mark, "MSH6_|SETD2_"),
                                                                "MMR",
                                                                ifelse(str_detect(dna_repair_mark, "XRCC4"),
                                                                       "DSBR",
                                                                       NA)))))
condition_pathway_pairs = metadata %>% 
  select(info1, info2) %>% 
  distinct %>% 
  filter(! str_detect(info2, "[C,c]ontrol|[O,o]ther")) %>% 
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
                                                                                 str_detect(info1, "Etoposide|Bleomycin|Camptothecin|Olaparib|Temozolomide|Melphalan|Cyclophosphamide|Mechlorethamine"),
                                                                               "DSBR",
                                                                               NA)))))) %>% 
  separate_rows(putative_repair_pathway_involved, sep = ",") %>% 
  filter(!is.na(putative_repair_pathway_involved)) %>% 
  left_join(metadata) %>% 
  filter(Sample %in% results_regressions$sample_id) %>% 
  arrange(putative_repair_pathway_involved) %>% 
  left_join(repair_mark_pathways) %>% 
  select(Sample, dna_repair_mark) %>% 
  rename("DNA repair\nactivity" = "dna_repair_mark") %>% 
  group_by(Sample) %>%
  mutate(Sample_possible_marks = n()) %>%
  ungroup()
condition_pathway_pairs = condition_pathway_pairs %>% 
  rowwise %>% 
  mutate(max_sample_mark_score = 1 / length(unique(condition_pathway_pairs$Sample)) / Sample_possible_marks) %>%
  select(-c(Sample_possible_marks))
condition_pathway_pairs = condition_pathway_pairs %>% 
  rowwise %>% 
  mutate(max_sample_mark_score = max_sample_mark_score / sum(condition_pathway_pairs$max_sample_mark_score))


## I actually run all k´s possible
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
    mutate(Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h))))) %>%
    # highlight top hits for each signature
    group_by(Signature) %>% 
    mutate(is.hit = ifelse(Exposure==max(Exposure), "Top hit", NA))
  
  # DNA repair mark weights in signatures
  weights = nmf_res$w %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column("dna_repair_mark") %>%
    pivot_longer(cols = contains("nmf"), names_to = "signature", values_to = "weight") %>% 
    extract(dna_repair_mark, into = c("dna_repair_mark", "submatrix"), "(.*)_([^_]+$)") %>% 
    mutate(dna_repair_mark = sub(".L", "", dna_repair_mark),
           dna_repair_mark = factor(dna_repair_mark, levels = unique(gsub(".L_...coeff", "", rownames(nmf_res$w)))),
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
    merge(condition_pathway_pairs) %>% 
    merge(exposures) %>% 
    select(Signature, Sample, Weight, Exposure, "DNA repair\nactivity", max_sample_mark_score) %>% 
    mutate(sample_mark_score = Weight * Exposure * max_sample_mark_score)

  write_tsv(good_model_score_table, "plots/good_model_score_table.tsv")
  
  good_model_score = sum(good_model_score_table$sample_mark_score)
  #####
  
  
  ### plotting

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
                  `colnames<-`("Signature"))}
  
  ## exposures (only samples with Exposures > exposure_limit%)
  pos = position_jitter(w = 0.25, h = 0, seed = 1)
  exposures_plot = ggplot(filtered_exposures, 
                          aes(x = Signature,
                              y = Exposure*100,
                              shape = dataset)) +
    scale_y_continuous(labels = function(x) sub("0+$", "", x)) +
    geom_point(aes(fill = info2),
               size = 4,
               position = pos) +
    scale_fill_manual(values = jet.colors(length(unique(exposures$info2)))) +
    scale_shape_manual(values = c(21, 24)) +
    guides(fill = guide_legend(override.aes = list(size=6, shape=21)),
           shape = guide_legend(override.aes = list(size=6))) +
    ggrepel::geom_text_repel(aes(label = paste(Sample, info1)),
                             size = 4,
                             force = 10,
                             position = pos,
                             max.overlaps = 1000000,
                             min.segment.length = 1) +
    facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
    theme_classic() +
    xlab("") +
    ggtitle(paste0("Model Score = ", round(good_model_score*100, 2), "%")) +
    ylab(paste0("% Exposure (>", exposure_limit*100, "%)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(4, "mm"),
          legend.title = element_blank(),
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
