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


#### Here we load all bin-wise mutcounts of all chromsomes CONCATENATED
# so we just sum up the mutcounts grouping by bin profile


# from command parameters
args = commandArgs(trailingOnly=TRUE)


dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[1]) %>%
  read_csv(comment = "#")


dnarep_mark_simulate = ifelse(interactive(),
                yes = "OGG1_GOx30_chipseq",
                no = args[2])

altered_level = ifelse(dnarep_mark_simulate == "AID_regions",
                       yes = "targets",
                       no = "high_activity")


mutfoldinc = ifelse(interactive(),
                yes = "2",
                no = args[3]) %>% as.numeric


## load ready_for_regression (ALL sep chromosomes) from previous process, and bind them
sim_pos_con = ifelse(interactive(),
                yes = lapply(list(c(Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr1_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr2_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr3_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr4_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr5_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr6_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr7_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr8_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr9_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr10_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr11_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr12_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr13_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr14_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr15_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr16_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr17_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr18_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr19_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr20_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr21_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chr22_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_sim_pos_con_chrX_mutfoldinc2_OGG1_GOx30_chipseq.tsv")[1]
                )),
                read_tsv),
                no = lapply(list(args[-(1:3)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)
gc()

simulated_mutcount_colname = names(sim_pos_con)[1]


trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")

sim_pos_con = sim_pos_con %>% 
  select(-simulated_mark) %>% 
  ## sum muts across all concatenated chromosomes bin-wise (i.e. grouped)
  group_by_at(vars(!matches(simulated_mutcount_colname))) %>% 
  summarise_at(simulated_mutcount_colname, ~sum(.)) %>%
  relocate(simulated_mutcount_colname) %>% 
  ungroup() %>% 
  # mb_domain and tri as ordered and unordered factors, respectively
  mutate(mb_domain = factor(mb_domain, ordered = T),
         tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
  arrange(tri, mb_domain) %>%
  as_tibble
gc()


### regression

# formula specifies which of the fold values that were used to simulate the mutcounts to use as dependent variable
formula = paste0(paste0("simulated_mutcount_", mutfoldinc, "fold ~ "),
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "mb_domain + ",
                 "tri + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear model for the negative binomial family
cat(sprintf('Running regression...\n'))
y = suppressWarnings(glm.nb(formula = formula, 
                            data = sim_pos_con))

## parse output
y_tidy = broom::tidy(y) %>%
  filter(term != "(Intercept)" & !str_detect(term, "mb_domain") & !str_detect(term, "^tri")) %>%
  # calc CI 95% from std error
  mutate("conf.low" = estimate - `std.error`*1.96,
         "conf.high" = estimate + `std.error`*1.96) %>% 
  select(c(term, estimate, contains("conf"))) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(sample_id = paste0(dnarep_mark_simulate, "-", altered_level, "__muts_", mutfoldinc, "-fold"),
         info1 = sample_id,
         info2 = sample_id,
         # theta value used, either because it was optimal or because it reached 10 iterations
         theta = y$theta) %>% 
  relocate(sample_id)
gc()

write_tsv(y_tidy, "simulated_positive_control.tsv")
