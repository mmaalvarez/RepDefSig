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


#### Here we load all bin-wise mutcounts of all chromosomes of the baseline sample CONCATENATED
# so we just sum up the mutcounts grouping by bin profile
# and create sim_pos_con
# and regression


# from command parameters
args = commandArgs(trailingOnly=TRUE)

# load utils.R (functions) -- only used here if trinucleotide matching/adj
if(interactive()){
  source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")
} else {
  source(args[1])
}

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[2]) %>%
  read_csv(comment = "#")


dnarep_mark_simulate = ifelse(interactive(),
                yes = "OGG1_GOx30_chipseq",
                no = args[3])

altered_level = ifelse(dnarep_mark_simulate %in% c("CTCF_cohesin","A3A_TpCpH_hairpins"),
                       yes = "targets",
                       no = "high_activity")

mutfoldinc = ifelse(interactive(),
                yes = "0.01",
                no = args[4]) %>% as.numeric

# load offset from 3rd process
offset = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1],
                no = args[5]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"
offset = offset %>% 
  mutate(mb_domain = factor(mb_domain))

## type of trinuc modification, if any
trinuc_mode = ifelse(interactive(),
                     yes = "adjustment",
                     no = args[6])


## load ready_for_regression (ALL sep chromosomes) (ACTUALLY BASELINE SAMPLE) from previous process, and bind them
baseline_sample = ifelse(interactive(),
                yes = lapply(list(c(Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr1_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr2_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr3_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr4_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr5_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr6_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr7_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr8_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr9_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr10_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr11_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr12_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr13_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr14_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr15_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr16_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr17_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr18_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr19_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr20_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr21_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chr22_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1],
                                    Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/baseline_sample_chrX_mutfoldinc0.01_OGG1_GOx30_chipseq.tsv")[1])),
                read_tsv),
                no = lapply(list(args[-(1:6)][str_detect(args[-(1:6)], 
                                                         paste0("_mutfoldinc", mutfoldinc, "_", dnarep_mark_simulate, ".tsv"))]), 
                            read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>% 
  ### ensure that we keep only mutfold and simulated feature specified in this nf channel
  filter(`mutfoldinc`==mutfoldinc & `simulated_mark`==dnarep_mark_simulate) %>% 
  select(-c(chr,mutfoldinc)) %>% 
  ### sum up the mutcounts across bound chromosomes, grouping by bin profile
  group_by(across(!contains("mean_mutcount") & !matches("chr"))) %>% 
  summarise_at(vars(contains("mean_mutcount")),
               ~ceiling(sum(.))) %>% 
  ungroup

gc()


#####

## creating sim_pos_con from baseline sample's merged chromosomes

## get total mut burden at background bin in baseline sample KEEP TRINUC SEP
bg_mut_burden = baseline_sample %>% 
  # get only mut counts in bg bin ("low" or "bgGenome")
  filter(get(dnarep_mark_simulate) %in% c("low", "bgGenome")) %>% 
  group_by(tri, mb_domain) %>% 
  # collapse bins
  summarise(mutfoldinc_mutcount = ceiling(sum(mean_mutcount) * mutfoldinc)) %>% 
  ungroup


### simulate pos control samples out of the baseline_sample

# store for later adding back
baseline_sample_low_bins = baseline_sample %>%
  filter(get(dnarep_mark_simulate) %in% c("low", "bgGenome")) %>% 
  rename("simulated_mutcount_x-fold" = mean_mutcount)

sim_pos_con = baseline_sample %>%
  # keep only high bins to speed up, then the unaltered low bins are re-added
  filter(get(dnarep_mark_simulate) %in% c("high", "CTCF_cohesin_peak", "hairpin_TpCpH")) %>% 
  #### increase mutations in bins that are i) "high" repair mark abundance -to approach the baseline-, or ii) "CTCF_cohesin_peak"|"hairpin_TpCpH" -to simulate an CTCF_cohesin_peak|A3A-SHM sample-
  ## mutfoldinc is multiplied to the baseline sample's genome-wide mutation burden per trinuc (e.g. if baseline sample's total # mutations at ATA in target bin (e.g. CTCF_cohesin_peak-targets) is just 2, and in bg bin (bgGenome) is 250, then multiply to the target bin a FRACTION of the total # ATA muts of the bg bin: e.g. a 0.01 --> 2 * 0.01 * 250
  # ceiling the mutfoldinc*mutburden ensures that it's not <1 in case mutburden was too low
  # mean_mutcount+1 ensures that bins with 0 muts get something
  merge(bg_mut_burden, all=T) %>% 
  mutate(`simulated_mutcount_x-fold` = ceiling((`mean_mutcount`+1) * mutfoldinc_mutcount)) %>%
  select(names(baseline_sample_low_bins)) %>% 
  bind_rows(baseline_sample_low_bins) %>% 
  # here we do leave the offset column
  merge(offset, all = T) %>%
  select(`simulated_mutcount_x-fold`,
         all_of(dnarep_marks$name),
         mb_domain,
         tri,
         trinuc32,
         freq_trinuc32) %>%
  mutate(simulated_mark = dnarep_mark_simulate)
gc()
########


## prepare tables for trinuc matching or adj, if selected
if(trinuc_mode %in% c("matching", "adjustment")){
  
  sim_pos_con_tmp = sim_pos_con %>% 
    select(-c(tri, simulated_mark)) %>% 
    unite("bin", !contains("mutcount") & !contains("trinuc32"))
  
  # table with total trinucleotide counts per bin
  total_trinuc_table = sim_pos_con_tmp %>% 
    select(-contains("mutcount")) %>% 
    distinct() %>% 
    pivot_wider(names_from = trinuc32, values_from = freq_trinuc32) %>% 
    column_to_rownames("bin")
  
  # table with MUTATED trinucleotide counts per bin
  mut_trinuc_table = sim_pos_con_tmp %>% 
    select(-freq_trinuc32) %>% 
    group_by(bin, trinuc32) %>% 
    summarise(across(contains("mutcount"), ~sum(.x), .names = "mutcount")) %>% 
    pivot_wider(names_from = trinuc32, values_from = mutcount) %>% 
    column_to_rownames("bin")
  
  rm(sim_pos_con_tmp) ; gc()
}

if(trinuc_mode == "matching"){
  
  matched_counts = trinuc_matching(total_trinuc_table, mut_trinuc_table,
                                   stoppingCriterion = 0.01,
                                   maxIter = 20000*length(sim_pos_con),
                                   n_finish_tokens = 1000,
                                   max_fraction_removed_muts = 0.25)
  
  matched_totalcounts = matched_counts$total_trinuc %>% 
    pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"), 
                 names_to = "trinuc32", values_to = "freq_trinuc32") %>% 
    mutate(bin = gsub("CTCF_cohesin_peak", "CTCFcohesinpeak", bin),
           bin = gsub("hairpin_TpCpH", "hairpinTpCpH", bin)) %>% 
    separate(bin, into = select(sim_pos_con, -c(simulated_mark, contains("mutcount"), tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("CTCFcohesinpeak", "CTCF_cohesin_peak", .)) %>% 
    mutate_if(is.character, ~gsub("hairpinTpCpH", "hairpin_TpCpH", .)) %>% 
    replace_na(list(freq_trinuc32 = 0)) %>% 
    # just in case there was some trinuc freq slightly below 0, make it 0
    mutate(freq_trinuc32 = ifelse(freq_trinuc32 <= -1,
                                  0,
                                  freq_trinuc32)) %>% 
    merge(distinct(select(offset, -c(tri,freq_trinuc32))), 
          all = T) %>% 
    replace_na(list(freq_trinuc32 = 0))
  
  reg_table = matched_counts$mut_trinuc %>% 
    pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"), 
                 names_to = "trinuc32", values_to = "mutcount") %>% 
    mutate(bin = gsub("CTCF_cohesin_peak", "CTCFcohesinpeak", bin),
           bin = gsub("hairpin_TpCpH", "hairpinTpCpH", bin)) %>% 
    separate(bin, into = select(sim_pos_con, -c(simulated_mark, contains("mutcount"), tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("CTCFcohesinpeak", "CTCF_cohesin_peak", .)) %>% 
    mutate_if(is.character, ~gsub("hairpinTpCpH", "hairpin_TpCpH", .)) %>% 
    merge(matched_totalcounts, all = T) %>%
    replace_na(list(mutcount = 0)) %>% 
    # convert offset to logarithmic
    mutate(log_freq_trinuc32 = log(freq_trinuc32 + 1)) %>%
    select(-c(trinuc32, freq_trinuc32)) %>%
    relocate(mutcount) %>%
    relocate(mb_domain, .before = "log_freq_trinuc32") %>%
    mutate(mb_domain = factor(mb_domain, ordered = T)) %>% 
    arrange(mb_domain) %>% 
    as_tibble
  
  formula = paste0("mutcount ~ ",
                   paste(dnarep_marks$name, collapse = " + "), " + ",
                   "mb_domain + ",
                   "offset(log_freq_trinuc32)")
  
} else if(trinuc_mode == "adjustment"){
  
  reg_table = trinuc_adjustment(total_trinuc_table, mut_trinuc_table) %>% 
    mutate(bin = gsub("CTCF_cohesin_peak", "CTCFcohesinpeak", bin),
           bin = gsub("hairpin_TpCpH", "hairpinTpCpH", bin)) %>% 
    separate(bin, into = select(sim_pos_con, -c(simulated_mark, contains("mutcount"), tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("CTCFcohesinpeak", "CTCF_cohesin_peak", .)) %>% 
    mutate_if(is.character, ~gsub("hairpinTpCpH", "hairpin_TpCpH", .)) %>%
    merge(distinct(select(offset, -c(tri, contains("trinuc32")))), 
          all = T) %>% 
    replace_na(list(mutcount = 0)) %>% 
    relocate(mutcount) %>%
    relocate(mb_domain, .after = last_col()) %>%
    mutate(mb_domain = factor(mb_domain, ordered = T)) %>% 
    arrange(mb_domain) %>% 
    as_tibble
  
  formula = paste0("mutcount ~ ",
                   paste(dnarep_marks$name, collapse = " + "), " + ",
                   "mb_domain")
  
} else { ### no trinuc modification
  
  trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")
  
  reg_table = sim_pos_con %>%
    select(-c(simulated_mark, freq_trinuc32)) %>% 
    merge(offset, all = T) %>% 
    rename("mutcount" = "simulated_mutcount_x-fold") %>%
    replace_na(list(mutcount = 0)) %>%
    # convert offset to logarithmic
    mutate(log_freq_trinuc32 = log(freq_trinuc32 + 1)) %>% 
    select(-c(trinuc32, freq_trinuc32)) %>% 
    relocate(mutcount) %>%
    relocate(mb_domain, .before = "tri") %>%
    # mb_domain and tri as ordered and unordered factors, respectively
    mutate(mb_domain = factor(mb_domain, ordered = T),
           tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
    arrange(tri, mb_domain) %>% 
    as_tibble
  
  formula = paste0("mutcount ~ ",
                   paste(dnarep_marks$name, collapse = " + "), " + ",
                   "mb_domain + ",
                   "tri + ",
                   "offset(log_freq_trinuc32)")
}
rm(sim_pos_con) ; gc()


### regression

if(sum(reg_table$mutcount) >= 1){ # only do regression if there are not only 0 muts
  
  # dna rep mark levels as factors
  reg_table = reg_table %>% 
    mutate_at(vars(contains(match = dnarep_marks$name)),
              ~if('CTCF_cohesin_peak' %in% unique(.)){
                factor(., ordered = F, levels = c('bgGenome', 'CTCF_cohesin_peak')) # x-axis:
                # #SNVs |
                #       | ------ <-- CTCFcohesinpeak-SHM in tumor (maybe not flat, but with less negative coeff.)
                #       | \
                #       |  \  <-- no CTCFcohesinpeak-SHM in tumor
                #       |___\_____ 
                #        bg  CTCF_cohesin_peaks
              }else if('hairpin_TpCpH' %in% unique(.)){
                factor(., ordered = F, levels = c('bgGenome', 'hairpin_TpCpH')) # x-axis:
                # #SNVs |
                #       | ------ <-- A3A expressed in tumor (maybe not flat, but with less negative coeff.)
                #       | \
                #       |  \  <-- A3A not expressed in tumor
                #       |___\_____ 
                #        bg  TpCpH_hairpins
              }else{
                factor(., ordered = F, levels = c('low', 'high'))}) # x-axis:
                # #SNVs |
                #       | ------ <-- BERdef tumor (maybe not flat, but with less negative coeff.)
                #       | \
                #       |  \  <-- BERwt tumor
                #       |___\_____ 
                #        low  high
                #        OGG1 OGG1

  # stick to a generalized linear model for the negative binomial family
  cat(sprintf('Running regression in trinucl. mode: "%s"...\n', trinuc_mode))
  y = suppressWarnings(glm.nb(formula = formula, 
                              data = reg_table))
  
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
  
  } else {
    # since there are no mutations, create an empty "y_tidy"
  
    cat(sprintf('WARNING: sample %s has 0 mutations: can not run regression...\n', paste0(dnarep_mark_simulate, "-", altered_level, "__muts_", mutfoldinc, "-fold")))
    
    col_names_orig = paste0(c(paste("estimate", dnarep_marks$name, sep = "_"), paste("conf.low", dnarep_marks$name, sep = "_"), paste("conf.high", dnarep_marks$name, sep = "_")), "high")
    col_names_mod = c("sample_id", 
                      gsub("CTCF_cohesinhigh", "CTCF_cohesinCTCF_cohesin_peak", col_names_orig),
                      "theta")
    col_names = gsub("A3A_TpCpH_hairpinshigh", "A3A_TpCpH_hairpinshairpin_TpCpH", col_names_mod)
    
    y_tidy = data.frame(matrix(ncol = length(col_names), nrow = 1)) %>% 
      `colnames<-`(col_names) %>% 
      mutate(sample_id = paste0(dnarep_mark_simulate, "-", altered_level, "__muts_", mutfoldinc, "-fold"),
             info1 = sample_id,
             info2 = sample_id) %>% 
      relocate(theta, .after = "info2")
}
gc()

write_tsv(y_tidy, "simulated_positive_control.tsv")
