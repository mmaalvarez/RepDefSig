library(tidyverse)
library(data.table)
library(rlang)
library(MASS)
library(lme4)
library(broom.mixed)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")


### append offset; regression

# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample = ifelse(interactive(),
                yes = "MSM0.103", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

metadata_sample = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv",
                         no = args[2]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info"))) %>% 
  filter(sample_id == sample)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[3]) %>%
  read_csv(comment = "#")


# load offset from 3rd process
offset = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1],
                no = args[4]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"


## load ready_for_regression (ALL collected chromosomes) from previous process
merged = ifelse(interactive(),
                yes = "ready_for_regression_chr21.tsv", #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr21.tsv")[1],
                no = args[5]) %>% 
  fread %>% 
  as_tibble
gc()


## add mutation counts per trinuc×mb_domain
trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")

sommut_tricount_dnarep_chromatin = merged %>%
  select(-c(mut_id, `.overlap`)) %>% 
  table %>%
  as.data.frame %>%
  rename("mutcount" = "Freq") %>% 
  ## add offset
  merge(offset, all = T) %>%
  replace_na(list(mutcount = 0)) %>% 
  relocate(mutcount) %>%
  relocate(mb_domain, .before = "log_freq_trinuc32") %>%
  relocate(tri, .after = "mb_domain") %>%
  # mb_domain and tri as ordered and unordered factors, respectively
  mutate(mb_domain = factor(mb_domain, ordered = T),
         tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
  arrange(tri, mb_domain) %>% 
  as_tibble

rm(merged) ; gc()


### regression

formula = paste0("mutcount ~ ",
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "(1 | mb_domain) + ",
                 "(1 | tri) + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear mixed-effects model for the negative binomial family
y = suppressWarnings(glmer.nb(formula = formula, 
                              data = sommut_tricount_dnarep_chromatin))

## parse output
y_tidy = broom.mixed::tidy(y, exponentiate = F, effects = "fixed") %>%
  # calc CI 95% from std error
  mutate("conf.low" = estimate - `std.error`*1.96,
         "conf.high" = estimate + `std.error`*1.96) %>% 
  select(c(term, estimate, contains("conf"))) %>% 
  filter(term != "(Intercept)") %>%
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(sample_id = sample,
         # theta value used, either because it was optimal or because it reached 10 iterations
         theta = lme4:::getNBdisp(y))
gc()

## append features' coefficients and pvalues to metadata_sample
results_real_sample = full_join(metadata_sample, y_tidy) %>%
  relocate(sample_id) %>% 
  relocate(info1, info2, .before = theta)
gc()

write_tsv(results_real_sample, "results_real_sample.tsv")
