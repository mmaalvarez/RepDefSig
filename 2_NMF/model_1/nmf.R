library(tidyverse)
library(RcppML)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")


# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../../1_parser_and_regressions/model_1/res/results.tsv") %>%
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
                           by = "sample_id") %>%
  # transpose
  t() %>% 
  `colnames<-`(.[1, ]) %>% 
  .[-1, ] %>%
  as_tibble %>% 
  mutate_all(as.numeric) %>%
  as.matrix

rm(results_regressions) ; rm(results_regressions_posmatrix) ; rm(results_regressions_negmatrix) ; gc()


## run NMF
nmf_res = nmf(coefficient_matrix,
              k = 35, # in the package paper (biorxiv.org/content/10.1101/2021.09.01.458620v1.full.pdf) this is the optimal number of signatures
              seed = 1)
# nmf_res$h --> Signature exposures in samples
# nmf_res$w --> DNA repair mark weights in signatures
