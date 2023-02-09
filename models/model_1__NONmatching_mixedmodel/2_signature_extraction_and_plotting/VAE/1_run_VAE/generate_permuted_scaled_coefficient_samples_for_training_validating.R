library(tidyverse)
library(data.table)
library(msm) # rtnorm
library(scales)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


# load original and positive control simulated coefficient tables
original_samples_table = read_tsv("../../../1_parser_and_regressions/res/results.tsv") %>% 
  select(sample_id, contains("estimate_"), contains("conf")) %>% 
  rename_all(~str_replace_all(., '.L', '')) %>% 
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat)

simulated_positive_controls_table = read_tsv("../../../1_parser_and_regressions/res/results_simulated_positive_controls.tsv") %>% 
  select(sample_id, contains("estimate_"), contains("conf")) %>% 
  rename_all(~str_replace_all(., '.L', '')) %>% 
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat)

# write the original (+ simulated pos controls) coefficients for the final processing with the trained VAE (to get signatures)
original_and_simposcon_scaled = original_samples_table %>% 
  bind_rows(simulated_positive_controls_table) %>% 
  select(!contains("conf")) %>% 
  group_by(mark) %>% 
  # scale the CI-permuted coefficients matrices to [-1, +1], because final decoding layer uses tanh as activation
  mutate(scaled_estimate = rescale(estimate,
                                   to = c(-1, 1))) %>% 
  ungroup %>% 
  select(-estimate) %>% 
  pivot_wider(names_from = mark, values_from = scaled_estimate)

# sample names must be in the first column, all other columns must be the numerical features:
#   sample_short	feature_1	feature_2	any_name ...
#     sample1	      -4	        3.4	      4
#     any_name	     4	        3.2	      0

write_tsv(original_and_simposcon_scaled, "VAE_input/original_and_simposcon_scaled.tsv")



##############################################################################

#### Run coefficient permutations (based on CI 95%) to generate samples for VAE training + validation

## Parameters and initializing of some objects

# number of permuted-samples tables (1 per epoch in VAE)
epochs = 200
# number of permuted samples per table
totalNumIters = 1000

coefficient_Resamp = list()
set.seed(1)


## Function to generate matrices resampling betas from their CI95% distributions (instead of UPmultinomial)
resample_from_CI_and_scale = function(original_samples_table){
  
  original_samples_table %>% 
    group_by(sample_id, mark) %>% 
    summarise(resampled_estimate = rtnorm(n = 1,
                                          mean = estimate,
                                          sd = 1,
                                          lower = conf.low,
                                          upper = conf.high)) %>% 
    # scale the CI-permuted coefficients matrices to [-1, +1] (e.g. python minmax scaler), because final decoding layer uses tanh as activation
    select(!contains("conf")) %>% 
    group_by(mark) %>% 
    mutate(scaled_estimate = rescale(resampled_estimate,
                                     to = c(-1, 1))) %>% 
    ungroup %>% 
    select(-resampled_estimate) %>% 
    pivot_wider(names_from = mark, values_from = scaled_estimate)
}

# run permutations
for(epoch in 1:epochs){
  
  for(nIter in 1:totalNumIters) {
    
    print(paste0("Generating permuted sample ", nIter, " for epoch ", epoch))
    
    # for each sample (row) resample coefficients from CI95% distrs, and scale
    coefficient_Resamp[[nIter]] = resample_from_CI_and_scale(original_samples_table) %>% 
      mutate("nIter" = nIter)
    
    gc()
  }
  
  # for a given epoch, bind all permuted samples as a VAE training+validation data input
  coefficient_Resamp = bind_rows(coefficient_Resamp) %>% 
    # shuffle all samples
    slice(sample(1:n()))
  
  write_tsv(coefficient_Resamp, 
            paste0("VAE_input/permuted_coefficients_", totalNumIters, "__epoch_", epoch, ".tsv"))
  coefficient_Resamp = list()
  gc()
}

