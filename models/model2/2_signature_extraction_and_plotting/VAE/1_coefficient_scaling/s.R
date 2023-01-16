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


######################
#### Generate matrices resampling coefficients from their CI95% distributions (instead of UPmultinomial)

coefficient_table = read_tsv("../../2_coefficient_permutations/original_data.tsv")

## write the original coefficients for the final processing with the trained VAE (to get signatures)
# scale the CI-permuted coefficients matrices to [-1, +1] (e.g. python minmax scaler), because final decoding layer uses tanh as activation
original_data_scaled = coefficient_table %>% 
  select(!contains("conf")) %>% 
  group_by(mark) %>% 
  mutate(scaled_estimate = rescale(estimate,
                                   to = c(-1, 1))) %>% 
  ungroup %>% 
  select(-estimate) %>% 
  pivot_wider(names_from = mark, values_from = scaled_estimate)

write_tsv(original_data_scaled, "original_data_scaled.tsv")


# load and scale resampled coefficients
totalNumIters = 1000
coefficient_Resamp_scaled = read_tsv(paste0("../../2_coefficient_permutations/permuted_coefficients_", totalNumIters, "iters.tsv")) %>% 
  group_by(mark, nIter) %>% 
  mutate(scaled_resampled_estimate = rescale(resampled_estimate,
                                             to = c(-1, 1))) %>% 
  ungroup %>% 
  select(-resampled_estimate) %>% 
  pivot_wider(names_from = mark, values_from = scaled_resampled_estimate) %>% 
  # reshuffle samples (i.e. row-wise) to maybe avoid biases in VAE training
  sample_frac(1L) %>% 
  select(-nIter)

write_tsv(coefficient_Resamp_scaled, paste0("VAE_input_", totalNumIters, "iters.tsv"))
          