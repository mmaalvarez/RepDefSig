library(tidyverse)
library(rlang)
library(MASS)
library(lme4)
library(broom)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")


### collect chroms, append offset; regression

# from command parameters
args = commandArgs(trailingOnly=TRUE)

metadata = ifelse(interactive(),
                  yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv",
                  no = args[1]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info")))

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[2]) %>%
  read_csv(comment = "#")

# load offset from 3rd process
offset = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1],
                no = args[3]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"


## load ready_for_regression (ALL sep chromosomes) from previous process, and bind them
merged = ifelse(interactive(),
                yes = lapply(list(c(Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr1.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr2.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr3.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr4.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr5.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr6.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr7.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr8.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr9.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr10.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr11.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr12.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr13.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr14.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr15.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr16.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr17.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr18.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr19.tsv")[1] #,
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr20.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr21.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr22.tsv")[1],
                                    #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chrX.tsv")[1]
                )),
                read_tsv),
                no = lapply(list(args[-(1:3)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)

sample = unique(merged$sample_id)

metadata_sample = metadata %>% 
  filter(sample_id == sample)

gc()


sommut_tricount_dnarep_chromatin = merged %>%
  select(-c(sample_id, mut_id, `.overlap`)) %>% 
  table %>%
  as.data.frame %>%
  rename("mutcount" = "Freq") %>% 
  ## add offset
  merge(offset, all = T) %>%
  replace_na(list(mutcount = 0)) %>% 
  relocate(mutcount) %>%
  relocate(mb_domain, .before = "log_freq_trinuc32") %>%
  # mb_domain as ordered factor
  mutate(mb_domain = factor(mb_domain, ordered = T)) %>% 
  arrange(mb_domain) %>% 
  as_tibble

rm(merged) ; gc()


### regression
formula = paste0("mutcount ~ ",
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "mb_domain + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear model for the negative binomial family
cat(sprintf('Running regression...\n'))
y = suppressWarnings(glm.nb(formula = formula, 
                            data = sommut_dnarep_chromatin))

## parse output
y_tidy = broom::tidy(y) %>%
  filter(term != "(Intercept)" & !str_detect(term, "mb_domain")) %>%
  # calc CI 95% from std error
  mutate("conf.low" = estimate - `std.error`*1.96,
         "conf.high" = estimate + `std.error`*1.96) %>% 
  select(c(term, estimate, contains("conf"))) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(sample_id = sample,
         # theta value used, either because it was optimal or because it reached 10 iterations
         theta = y$theta)
gc()

## append features' coefficients and pvalues to metadata_sample
results_sample = full_join(metadata_sample, y_tidy) %>%
  relocate(sample_id) %>% 
  relocate(info1, info2, .before = theta)
gc()

write_tsv(results_sample, "results_sample.tsv")
