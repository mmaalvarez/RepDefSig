library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
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


# - For each DNA repair pathway (e.g. MMR):
# 	- 1st) calculate the mutation burden baseline (=="baseline sample"): mean # mutations (treating each of the 96 trinuc types separately, to keep relative %s thereof) across all control samples, WITHIN each of the 2 genome BINS (i.e. high/low MMR repair abundance)
# 	- 2nd) each simulated sample will be the "baseline sample" + an increased # mutations (×2, ×3, ×4…), but ONLY in the genome bin with 'high' activity of a mark for that repair (e.g. high MSH6)
# 		- again, each trinuc type treated independently to keep the relative trinuc %s
# - All the simulated samples will be then run together with real ones in the NMF/VAE
# - The simulated samples should have high exposures of the signatures with high weights for the expected mark (e.g. for "MMRdef simulated samples", high "MSH6-weight signature exposure")
# 	- the more sensitive the model, low mut. burden-increase (e.g. ×2?) simulated samples should already show high exposures


# from command parameters
args = commandArgs(trailingOnly=TRUE)

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_",
                                no = args[1]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

# controls for creating the 'baseline sample', for then adding to it extra mutations to create simulated positive controls
control_samples = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv",
                         no = args[2]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  map_df(~read_tsv(.x)) %>% 
  filter(str_detect(info2, "[C,c]ontrol")) %>% 
  pull(sample_id)

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
                yes = "../work/6d/8006f697193d072f9c3da5e6fe4fb9/offset.tsv",
                no = args[5]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"


## load map_features (all chromosomes) from 2nd process
dfleft = ifelse(interactive(),
                yes = "../res/map_features.tsv",
                no = paste0(args[6], "/res/map_features.tsv")) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()


## load mutation fold increases to apply on baseline sample
mutation_foldinc = ifelse(interactive(),
                          yes = "2,4,8,16", # i.e. mut burden at "high" bins multiplied by ×2, ×4...
                          no = args[7]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% as.numeric


## load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c("../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_OGG1_GOx30_chipseq.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_MSH6_control.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_OGG1_GOx60_chipseq.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_SETD2_control.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_UV_XRseq_NHF1_CPD_1h.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv",
                                           "../work/6d/8006f697193d072f9c3da5e6fe4fb9/median_score_XRCC4.tsv")), 
                                    read_tsv),
                       no = lapply(list(args[-(1:7)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)



## load control samples' somatic mutations

for(control_sample in control_samples){
  
  # find path+sample name that exists
  existing_file = c()
  for(file_exists in paste0(path_somatic_variation, control_sample, ".csv")){
    if(file.exists(file_exists)){
      existing_file = c(existing_file, file_exists)
    }
  }
  
  # raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
  if(length(existing_file) != 1){
    stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
  }
  
  ## load sample
  dfright = read_csv(existing_file) %>%
    select(chr, start, end, tri) %>% 
    # add +-50bp buffer around each SNV (i.e. 0.1kb windows), to ensure most of them overlap with reptime and DNA mark genomic ranges (I still lose quite a lot)
    mutate(start = start - 50,
           end = end + 50) %>%
    rename("chrom" = "chr") %>%
    mutate(mut_id = paste0("mut_", row_number()))
  gc()

  ## map chromatin features
  merged = bed_intersect(dfleft, dfright, suffix = c("_dfleft",
                                                     "_dfright")) %>%
    select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
    ## weighted averages of chromatin feature bins and dna repair abundances
    mutate(`.overlap` = `.overlap` + 1) %>%
    group_by(mut_id_dfright, tri_dfright) %>%
    # first get total bp of genome chunks with RT + dna marks info with which the SNV window overlap (max. 101)
    summarise(total_overlaps_length = sum(`.overlap`), across()) %>%
    # now calculate weighted (based on chunk length by total overlapping length of the window) RT and dna marks scores
    mutate_at(vars(!contains("mut_id_dfright") & !contains("tri_dfright") & !contains("total_overlaps_length") & !contains(".overlap")),
              ~. * `.overlap` / total_overlaps_length) %>%
    # ... and sum them across each SNV (mut_id)
    summarise_at(vars(!contains("mut_id_dfright") & !contains("total_overlaps_length") & !contains("tri_dfright") & !contains(".overlap")),
                 ~sum(.)) %>%
    ungroup %>% 
    # round the reptime bins to 0 decimal
    mutate_at(vars(contains(chromatin_features$name)),
              ~round(., 0)) %>% 
    # remove the "dleft" and "dright" parts of the column names
    rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
    # combine chromatin features (although there should typically be only RepliSeq)
    unite("mb_domain", contains(chromatin_features$name), sep = "_") %>% 
    # binarize weighted average DNA repair value by being lower or larger than the across-genome median
    rowwise %>% 
    mutate_at(vars(contains(match = dnarep_marks$name)),
              function(x){var_name = rlang::as_label(substitute(x))
              ifelse(x <= filter(median_scores, dnarep_mark == var_name)$median_score,
                     "low",
                     "high")}) %>%
    ### mut count
    select(-mut_id) %>% 
    table %>%
    as.data.frame %>%
    rename("mutcount" = "Freq") %>% 
    # add offset to have 0 muts in the non-existing combinations
    merge(offset, all = T) %>%
    replace_na(list(mutcount = 0)) %>% 
    relocate(mutcount) %>% 
    select(-log_freq_trinuc32) %>% 
    mutate("control_sample" = control_sample)
  
  ### append to control_samples_table if this exists, otherwise create it
  if(exists("control_samples_table")){
    control_samples_table = bind_rows(control_samples_table, merged)
  }else{
    control_samples_table = merged
  }
  rm(merged)
  gc()
}


### create baseline sample

# first exclude control outliers (too many mutations, as maybe there are cryptic mutagens)
mut_burden_controls = control_samples_table %>% 
  group_by(control_sample) %>% 
  summarise(mutburden = sum(mutcount)) %>% 
  pull(mutburden) %>% 
  quantile
# assert that there are no "extreme" outliers
#(mut_burden_controls[["75%"]] + 1.5 * (mut_burden_controls[["75%"]] - mut_burden_controls[["25%"]]))  >  max(mut_burden_controls)

control_samples_no_outliers = control_samples_table %>% 
  group_by(control_sample) %>% 
  summarise(mutburden = sum(mutcount)) %>% 
  # use as cutoff for removing treated samples the mean of mut_burden_controls
  filter(mutburden <= mean(mut_burden_controls)) %>% 
  pull(control_sample)

baseline_sample = control_samples_table %>% 
  filter(control_sample %in% control_samples_no_outliers) %>% 
  select(-control_sample) %>% 
  group_by(across(!contains("mutcount"))) %>% 
  summarise(mean_mutcount = mean(mutcount)) %>% 
  ungroup


### simulate pos control samples out of the baseline_sample
simulations_by_mark = list()
trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")

for(dnarep_mark in dnarep_marks$name){
  
    sim_pos_con = baseline_sample %>% 
      select(tri, mb_domain, dnarep_mark, mean_mutcount) %>% 
      group_by(tri, mb_domain, !!sym(dnarep_mark)) %>% 
      # collapse bins
      summarise(sum_mean_mutcount = sum(mean_mutcount))
    
    for(mutfoldinc in mutation_foldinc){
      
      sim_pos_con = sim_pos_con %>% 
        # add mutations in "high"
        mutate("simulated_mutcount_{mutfoldinc}fold" := ifelse(get(dnarep_mark) == "high",
                                                               yes = ceiling(sum_mean_mutcount * mutfoldinc),
                                                               no = ceiling(sum_mean_mutcount)))
    }
    
    # parse output
    sim_pos_con = sim_pos_con %>% 
      merge(offset, all = T) %>% 
      select(starts_with("simulated_mutcount"),
             dnarep_marks$name,
             mb_domain,
             tri,
             log_freq_trinuc32) %>% 
      ungroup %>% 
      # dna rep mark levels as ordered factors
      mutate_at(vars(contains(match = dnarep_marks$name)),
                ~ factor(., ordered = T, levels = c('low', 'high'))) %>% 
      # mb_domain and tri as ordered and unordered factors, respectively
      mutate(mb_domain = factor(mb_domain, ordered = T),
             tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
      arrange(tri, mb_domain) %>% 
      as_tibble
    
    # append to list
    simulations_by_mark[[dnarep_mark]] = sim_pos_con
    gc()
}


### regressions

for(dnarep_mark in dnarep_marks$name){
  
  for(mutfoldinc in mutation_foldinc){
    
    # formula specifies which of the fold values that were used to simulate the mutcounts to use as dependent variable
    formula = paste0(paste0("simulated_mutcount_", mutfoldinc, "fold ~ "),
                     paste(dnarep_marks$name, collapse = " + "), " + ",
                     "(1 | mb_domain) + ",
                     "(1 | tri) + ",
                     "offset(log_freq_trinuc32)")
    
    # first try a generalized linear mixed-effects model for the negative binomial family
    y = tryCatch(glmer.nb(formula = formula, 
                          # specify a dnarep_mark
                          data = simulations_by_mark[[dnarep_mark]]),
                 # if there is a failure to converge to theta, run Poisson
                 warning = function(w) glmer(formula = formula, 
                                             data = simulations_by_mark[[dnarep_mark]], 
                                             family = poisson),
                 error = function(e) glmer(formula = formula, 
                                           data = simulations_by_mark[[dnarep_mark]], 
                                           family = poisson))
    # track whether NB or Poisson
    nb_or_pois = family(y)[[1]]
    
    ## parse output
    y_tidy = broom.mixed::tidy(y,
                               exponentiate = F, effects = "fixed") %>%
      # calc CI 95% from std error
      mutate("conf.low" = estimate - `std.error`*1.96,
             "conf.high" = estimate + `std.error`*1.96) %>% 
      select(c(term, estimate, contains("conf"))) %>% 
      filter(term != "(Intercept)") %>%
      pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
      mutate(sample_id = paste0(dnarep_mark, "-high_muts_", mutfoldinc, "-fold"),
             glm = ifelse(str_detect(nb_or_pois, "Negative Binomial"),
                          "Negative Binomial",
                          ifelse(str_detect(nb_or_pois, "[P,p]oisson"),
                                 "Poisson",
                                 stop(paste0("ERROR! GLM family used was not 'Negative Binomial' nor '[P,p]oisson'; instead, it was '", nb_or_pois, "'"))))) %>% 
      relocate(sample_id)
    
    # append to final table if this exists, otherwise create it
    if(exists("results_simulated_positive_controls")){
      results_simulated_positive_controls = bind_rows(results_simulated_positive_controls, y_tidy)
    }else{
      results_simulated_positive_controls = y_tidy
    }
    gc()
  }
}

write_tsv(results_simulated_positive_controls, "results_simulated_positive_controls.tsv")
