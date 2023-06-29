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

sample = ifelse(interactive(),
                yes = "AY9808_REP1", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

metadata = ifelse(interactive(),
                  yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv",
                  no = args[2]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(`sample_id`, starts_with("info")))

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[3]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[4]) %>%
  read_csv(comment = "#")

# load offset from 3rd process
offset = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1],
                no = args[5]) %>% 
  read_tsv
# rename the chromatin environment column (if there is any, typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
if(length(chromatin_features) == 1){
  colnames(offset)[1] = "mb_domain"
  offset = offset %>% 
    mutate(mb_domain = factor(mb_domain, ordered = T))
}

## load ready_for_regression (ALL sep chromosomes) from previous process, and bind them
all_chr_sample_of_interest = args[-(1:5)] %>% 
  data.frame() %>% `colnames<-`("file_path_chr_sample_of_interest") %>% 
  filter(str_detect(file_path_chr_sample_of_interest, sample)) %>% 
  pull(file_path_chr_sample_of_interest)

merged = ifelse(interactive(),
                yes = lapply(list(c(# Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr1.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr2.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr3.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr4.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr5.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr6.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr7.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr8.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr9.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr10.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr11.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr12.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr13.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr14.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr15.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr16.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr17.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr18.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr19.tsv"))[1],
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr20.tsv"))[1],
                                    Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr21.tsv"))[1],
                                    Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chr22.tsv"))[1]#,
                                    # Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/ready_for_regression_", sample, "_chrX.tsv"))[1]
                                    )),
                             read_tsv),
                no = lapply(list(all_chr_sample_of_interest), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)
gc()

metadata_sample = metadata %>% 
  filter(sample_id == sample)


trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")


reg_table = select(merged, -c(sample_id, mut_id, `.overlap`)) %>% 
  table %>%
  as.data.frame %>%
  rename("mutcount" = "Freq") %>% 
  ## add offset
  merge(offset, all = T) %>%
  replace_na(list(mutcount = 0)) %>% 
  relocate(tri, .before = "trinuc32") %>%
  select(-c(freq_trinuc32)) %>% 
  merge(offset, all = T) %>%
  replace_na(list(mutcount = 0)) %>%
  # convert offset to logarithmic
  mutate(log_freq_trinuc32 = log(freq_trinuc32 + 1)) %>% 
  select(-c(trinuc32, freq_trinuc32)) %>% 
  relocate(mutcount) %>%
  # tri as ordered factor
  mutate(tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
  arrange(tri) %>% 
  as_tibble

reg_variables = case_when(length(chromatin_features$name) == 1 ~ paste0(paste(dnarep_marks$name, collapse = " + "), " + mb_domain + "),
                          length(chromatin_features$name) == 0 ~ paste0(paste(dnarep_marks$name, collapse = " + "), " + "))
formula = paste0("mutcount ~ ",
                 reg_variables,
                 "tri + ",
                 "offset(log_freq_trinuc32)")

rm(merged) ; gc()


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
              }else if("low" %in% unique(.)  &  "high" %in% unique(.)){
                factor(., ordered = F, levels = c('low', 'high')) # x-axis:
                # #SNVs |
                #       | ------ <-- BERdef tumor (maybe not flat, but with less negative coeff.)
                #       | \
                #       |  \  <-- BERwt tumor
                #       |___\_____ 
                #        low  high
                #        OGG1 OGG1
              }else{ ## in case levels are not "low" and "high", just let them be as they are, probably alphanumerically
                factor(., ordered = F) # x-axis:
              })
  
  # generalized linear model for the negative binomial family
  cat(sprintf('Running regression...\n'))
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
  
  col_names_orig = paste0(c(paste("estimate", dnarep_marks$name, sep = "_"), paste("conf.low", dnarep_marks$name, sep = "_"), paste("conf.high", dnarep_marks$name, sep = "_")), "high")
  col_names_mod = c("sample_id", 
                    gsub("CTCF_cohesinhigh", "CTCF_cohesinCTCF_cohesin_peak", col_names_orig),
                    "theta")
  col_names = gsub("A3A_TpCpH_hairpinshigh", "A3A_TpCpH_hairpinshairpin_TpCpH", col_names_mod)
  
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
