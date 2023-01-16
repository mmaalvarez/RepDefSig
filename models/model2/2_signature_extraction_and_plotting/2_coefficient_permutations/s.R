library(tidyverse)
library(data.table)
library(parallel)
library(msm) # rtnorm
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


sample_ids = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
               "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  map_df(~read_tsv(.x)) %>% 
  pull(sample_id)

results_regressions = read_tsv("../../1_parser_and_regressions/res/results.tsv") %>%
  filter(sample_id %in% sample_ids) %>% 
  select(sample_id, contains("estimate_"), contains("conf"))

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

write_tsv(coefficient_table, "original_data.tsv")



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

# bind all permutated tables as an input for the NMF and VAE (the latter after -1,1 scaling)
coefficient_Resamp = bind_rows(coefficient_Resamp)

write_tsv(coefficient_Resamp, paste0("permuted_coefficients_", totalNumIters, "iters.tsv"))
