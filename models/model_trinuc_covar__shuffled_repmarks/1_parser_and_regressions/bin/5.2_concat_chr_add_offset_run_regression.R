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

# load utils.R (functions) -- only used here if trinucleotide matching/adj
if(interactive()){
  source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")
} else {
  source(args[1])
}

sample = ifelse(interactive(),
                yes = "b668939e-7d77-504c-abfa-7b4982106ab9", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[2])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

metadata = ifelse(interactive(),
                  yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv",
                  no = args[3]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info")))

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[4]) %>%
  read_csv(comment = "#")

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


## load ready_for_regression (ALL sep chromosomes) from previous process, and bind them
merged = ifelse(interactive(),
                yes = lapply(list(c(#Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr1.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr2.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr3.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr4.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr5.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr6.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr7.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr8.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr9.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr10.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr11.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr12.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr13.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr14.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr15.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr16.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr17.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr18.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr19.tsv"),
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr20.tsv"),
                  Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr21.tsv"),
                  Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chr22.tsv") #,
                  #Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_chrX.tsv")
                )),
                read_tsv),
                no = lapply(list(args[-(1:6)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>% 
  ## keep only the sample of this channel
  filter(`sample_id` == sample)
  
gc()

metadata_sample = metadata %>% 
  filter(sample_id == sample)

merged = select(merged, -c(sample_id, mut_id, `.overlap`)) %>% 
  table %>%
  as.data.frame %>%
  rename("mutcount" = "Freq") %>% 
  ## add offset
  merge(offset, all = T) %>%
  replace_na(list(mutcount = 0)) %>% 
  relocate(tri, .before = "trinuc32")
gc()


## prepare tables for trinuc matching or adj, if selected
if(trinuc_mode %in% c("matching", "adjustment")){
  
  merged_tmp = merged %>% 
    select(-tri) %>% 
    unite("bin", !contains("mutcount") & !contains("trinuc32"))
  
  # table with total trinucleotide counts per bin
  total_trinuc_table = merged_tmp %>% 
    select(-mutcount) %>% 
    distinct() %>% 
    pivot_wider(names_from = trinuc32, values_from = freq_trinuc32) %>% 
    column_to_rownames("bin") #%>% 
    ## REMOVE rows with all counts 0
    #filter_all(., any_vars(. != 0))
  
  # table with MUTATED trinucleotide counts per bin
  mut_trinuc_table = merged_tmp %>% 
    select(-freq_trinuc32) %>% 
    group_by(bin, trinuc32) %>% 
    summarise(mutcount = sum(mutcount)) %>% 
    pivot_wider(names_from = trinuc32, values_from = mutcount) %>% 
    column_to_rownames("bin") #%>% 
    #filter_all(., any_vars(. != 0))
  
  rm(merged_tmp) ; gc()
}

if(trinuc_mode == "matching"){
  
  matched_counts = trinuc_matching(total_trinuc_table, mut_trinuc_table,
                                   stoppingCriterion = 0.01,
                                   maxIter = 20000*length(merged),
                                   n_finish_tokens = 1000,
                                   max_fraction_removed_muts = 0.25)
  
  matched_totalcounts = matched_counts$total_trinuc %>% 
    pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"), 
                 names_to = "trinuc32", values_to = "freq_trinuc32") %>% 
    mutate(bin = gsub("AID_", "AID", bin)) %>% 
    separate(bin, into = select(merged, -c(mutcount, tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("AIDtarget", "AID_target", .)) %>% 
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
    mutate(bin = gsub("AID_", "AID", bin)) %>% 
    separate(bin, into = select(merged, -c(mutcount, tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("AIDtarget", "AID_target", .)) %>% 
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
    mutate(bin = gsub("AID_", "AID", bin)) %>% 
    separate(bin, into = select(merged, -c(mutcount, tri, contains("trinuc32"))) %>% names) %>% 
    mutate_if(is.character, ~gsub("AIDtarget", "AID_target", .)) %>%
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
  
  reg_table = merged %>%
    select(-c(freq_trinuc32)) %>% 
    merge(offset, all = T) %>%
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
rm(merged) ; gc()


### regression

if(sum(reg_table$mutcount) >= 1){ # only do regression if there are not only 0 muts
  
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
    mutate(sample_id = sample,
           # theta value used, either because it was optimal or because it reached 10 iterations
           theta = y$theta)
  
} else {
  # since there are no mutations, create an empty "y_tidy"
  
  cat(sprintf('WARNING: sample %s has 0 mutations: can not run regression...\n', sample))
  
  col_names = c(gsub("AID_regionslow", "AID_regionsbgGenome", paste0(c(paste("estimate", dnarep_marks$name, sep = "_"), paste("conf.low", dnarep_marks$name, sep = "_"), paste("conf.high", dnarep_marks$name, sep = "_")), "low")), "sample_id", "theta")
  
  y_tidy = data.frame(matrix(ncol = length(col_names), nrow = 1)) %>% 
    `colnames<-`(col_names) %>% 
    mutate(sample_id = sample)
}
gc()

## append features' coefficients and pvalues to metadata_sample
results_real_sample = full_join(metadata_sample, y_tidy) %>%
  relocate(sample_id) %>% 
  relocate(info1, info2, .before = theta)
gc()

write_tsv(results_real_sample, "results_real_sample.tsv")
