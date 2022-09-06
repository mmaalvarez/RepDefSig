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



##### parse input

metadata = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv") %>% 
  select(sample_id, source, tissue, OriginalType, hr_status, MSI_status, smoking_history, treatment_platinum, treatment_5FU, tumorPurity, gender) %>% 
  rename("Sample" = "sample_id")

# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../../1_parser_and_regressions/model_1/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_")) %>%  #  source, MSI_status, hr_status, smoking_history, treatment_platinum, treatment_5FU, tissue, OriginalType, 
  arrange(sample_id)

# a) keep positive coefficients, and convert negative to zero 
results_regressions_posmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))

# b) convert positive coefficients to zero, and convert negative to positive
results_regressions_negmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))

# merge them into the NMF input
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



#### evaluate best combination of n variables and k signatures

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
coefficient_matrix_Resamp = list()
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
  nmfOutputByIter = parLapply(cl = cl,
                              X = coefficient_matrix_Resamp,
                              fun = function(x, nFact) {
                                x = x[, colSums(x) > 0]
                                nmf(x, rank = nFact, maxIter = 10000, seed = attr(x, "nIter"))},
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
sigClusterScores = optimizeClustersPar(nmfHmatAllByFact, maxNumClusters = maxK)

## Visualize the clustering quality scores in heatmaps (facetted by the 'smooth' parameter)
heatmap_clustering = sigClusterScores[k <= maxK] %>% melt(id.vars = c("nFact", "k"),
                                                          measure.vars = "avgClScore" # "minClScore" or  "avgClScore" or "secondWorstClScore"
) %>%
  ggplot(aes(nFact, k)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.4)
ggsave("NMF_heatmap_clustering.pdf",
       plot = heatmap_clustering,
       device = "pdf",
       width = 12.5,
       height = 6,
       dpi = 600)

## FOR 6 features:
optimal_k = 6



###### now run the final NMF using the optimal k and n features

# regenerate the coeff matrix, without mult. by 100, and transposing, for RcppML::nmf()
coefficient_matrix_RcppML = merge(results_regressions_posmatrix,
                                  results_regressions_negmatrix,
                                  by = "sample_id",
                                  suffixes = c("_poscoeff", "_negcoeff")) %>%
  # transpose
  t() %>% 
  `colnames<-`(.[1, ]) %>% 
  .[-1, ] %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column("dna_repair_mark") %>% 
  mutate(dna_repair_mark = gsub("estimate_", "", dna_repair_mark)) %>% 
  column_to_rownames("dna_repair_mark") %>% 
  mutate_all(as.numeric) %>%
  as.matrix

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
  # add info on MSI, HR, smoking, and treatments
  left_join(metadata) %>% 
  mutate(Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h)))),
         MSI_parsed = ifelse(MSI_status %in% c("HYPER", "MSI", "ERCC2mut"), "MSI", NA),
         hr_parsed = ifelse(hr_status %in% c("HR_deficient"), "HRdef", NA),
         smoking_history_parsed = ifelse(smoking_history %in% c("Current", "Former"), "Smoker", NA),
         treatment_platinum_parsed = ifelse(treatment_platinum == "TRUE", "Platinum", NA),
         treatment_5FU_parsed = ifelse(treatment_5FU == "TRUE", "5FU", NA)) %>%
  unite(col = "Metadata", MSI_parsed, hr_parsed, smoking_history_parsed,treatment_platinum_parsed, treatment_5FU_parsed, na.rm = T, sep = " & ") %>% 
  mutate(Metadata = gsub("^$", "NOTA/NA", Metadata)) %>% 
  rename("Database" = "source") %>% 
  # highlight top hits for each signature
  group_by(Signature) %>% 
  mutate(is.hit = ifelse(Exposure==max(Exposure), "hit", NA),
         has.Metadata = ifelse(Metadata != "NOTA/NA", "yes", "no"))
# "NOTA/NA" to the end
exposures$Metadata = factor(exposures$Metadata, levels = c(unique(exposures$Metadata)[unique(exposures$Metadata) != "NOTA/NA"],
                                                           "NOTA/NA"))
write_tsv(exposures,
          "NMF_exposures.tsv")

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
  ## keep larger of the 2 values for the 2 nmf submatrices (1 based on pos coefficients submatrix, and another on the absolute neg coefficients submatrix)
  summarise(Weight = max(weight)) %>% 
  ungroup %>% 
  rename("DNA repair activity" = "dna_repair_mark", 
         "Signature" = "signature") %>% 
  relocate(Signature)
write_tsv(weights,
          "NMF_weights.tsv")


### plotting

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## exposures (only Exposures > 0.001)
exposures_plot = ggplot(filter(exposures, Exposure > 0.001), 
                        aes(x = Signature,
                            # convert exposures to %
                            y = Exposure*100,
                            group = Database)) +
  scale_y_log10(labels = function(x) sub("0+$", "", x)) +
  # all points, no colors
  geom_point(aes(shape = Database),
             position = position_dodge(width = 1),
             alpha = 0.5,
             size = 4) +
  # colored those with metadata info
  geom_point(aes(shape = Database,
                 fill = Metadata,
                 alpha = has.Metadata),
             position = position_dodge(width = 1),
             size = 4) +
  scale_shape_manual(values = c(21,23,24,22,20,25)) + #"\u2716"
  scale_fill_manual(values = c(jet.colors(length(unique(filter(exposures, Exposure > 0.001 & Metadata!="NOTA/NA")$Metadata))), "white")) + # Metadata==NOTA/NA is assigned white
  scale_alpha_manual(values = c(0, 0.8), guide = 'none') +
  guides(shape = guide_legend(override.aes = list(size=6)),
         fill = guide_legend(override.aes = list(size=6, shape=21))) +
  # label sample names to top hits
  ggrepel::geom_text_repel(data = filter(exposures, is.hit == "hit"),
                           aes(label = Sample),
                           size = 3,
                           min.segment.length = 10000) +
  facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
  theme_classic() +
  xlab("") +
  ylab("% Exposure (>0.1% ; log10 scale)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1, "mm"))

## weights
weights_plot = ggplot(weights %>%
                        mutate(`DNA repair activity` = gsub("_2strands", "", `DNA repair activity`),
                               Signature = factor(gsub("nmf", "", Signature), levels = seq(1:length(levels(weights$Signature))))), 
                      aes(x = Signature,
                          y = Weight)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) sub("0+$", "", x)) +
  geom_col(aes(fill = `DNA repair activity`)) +
  scale_fill_manual(values = jet.colors(length(levels(weights$`DNA repair activity`)))) +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_grid(cols = vars(Signature), scales = "free", space = "free") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1, "mm"))

combined_plots = cowplot::plot_grid(NULL,
                                    exposures_plot,
                                    NULL,
                                    weights_plot,
                                    nrow = 4,
                                    rel_heights = c(0.02, 1,-0.05,1))
ggsave("NMF_plot.pdf",
       plot = combined_plots,
       device = "pdf",
       width = 22.5,
       height = 12,
       dpi = 600,
       bg = "white")
