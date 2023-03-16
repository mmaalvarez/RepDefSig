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


sim_pos_con = sim_pos_con %>% 
  select(-simulated_mark) %>% 
  ## sum muts across all concatenated chromosomes bin-wise (i.e. grouped)
  group_by_at(vars(!matches(simulated_mutcount_colname))) %>% 
  summarise_at(simulated_mutcount_colname, ~sum(.)) %>%
  relocate(simulated_mutcount_colname) %>% 
  ungroup() %>% 
  # mb_domain as factor
  mutate(mb_domain = factor(mb_domain, ordered = T)) %>% 
  arrange(mb_domain) %>%
  as_tibble
gc()


### regression


# formula specifies which of the fold values that were used to simulate the mutcounts to use as dependent variable
formula = paste0(paste0("simulated_mutcount_", mutfoldinc, "fold ~ "),
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "mb_domain + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear model for the negative binomial family
cat(sprintf('Running regression...\n'))
y = suppressWarnings(glm.nb(formula = formula, 
                            data = sim_pos_con))

## parse output
y_tidy = broom::tidy(y) %>%
  filter(term != "(Intercept)" & !str_detect(term, "mb_domain")) %>%
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
