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


path_somatic_variation = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_" %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[4]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[5]) %>%
  read_csv(comment = "#")

offset = ifelse(interactive(),
                yes = "../work/b0/0125c17277870afcb86aaab4c93782/offset.tsv",
                no = args[6]) %>% 
  read_tsv
colnames(offset)[1] = "mb_domain"

dfleft = ifelse(interactive(),
                yes = "../res/map_features.tsv",
                no = paste0(args[7], "/res/map_features.tsv")) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")

median_scores = ifelse(interactive(),
                       yes = lapply(list(c("../work/26/08b67aa330177901dc72e70df9cc6b/median_score_OGG1_GOx30_chipseq.tsv",
                                           "../work/02/bd5527ea98ae2fbcb50d1c7947a0fb/median_score_OGG1_GOx60_chipseq.tsv",
                                           "../work/fb/0c57088546769b3d33c5ce901bdda8/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv",
                                           "../work/f2/acd0d55a2fb0d7a07260a120f0d01e/median_score_UV_XRseq_NHF1_CPD_1h.tsv",
                                           "../work/d2/65d1e94a5015f7106d266899bae572/median_score_XRCC4.tsv",
                                           "../work/58/2192db8e29c53114153463b867ce68/median_score_SETD2_control.tsv",
                                           "../work/73/4d3c3f1693923f3808026f1194d58f/median_score_MSH6_control.tsv",
                                           "../work/35/ba92f745984eec59ec4ca62cd3d2df/median_score_TP53_dauno_MOLM13.tsv",
                                           "../work/db/db4afaf45c9ec6f18e486dc5dfcc87/median_score_TP53_dauno_K562.tsv",
                                           "../work/fb/8caddf7603a8d872bd5df497010bea/median_score_AID_regions.tsv")), 
                                    read_tsv),
                       no = lapply(list(args[-(1:8)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)

trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")


failed_samples = read_tsv("../failed", col_names = F)$X1

list_failed_samples = list()
gc()

for(sample in failed_samples){
  
  existing_file = c()
  for(file_exists in paste0(path_somatic_variation, sample, ".csv")){
    if(file.exists(file_exists)){
      existing_file = c(existing_file, file_exists)
    }
  }
  if(length(existing_file) != 1){
    stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
  }
  
  dfright = read_csv(existing_file) %>%
    select(chr, start, end, tri) %>% 
    rename("chrom" = "chr") %>%
    mutate(mut_id = paste0("mut_", row_number()))
  
  merged = dfright %>%
    bed_intersect(dfleft, suffix = c("_dfright", "_dfleft")) %>%
    select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
    # remove the "dleft" and "dright" parts of the column names
    rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
    # combine chromatin features (although there should typically be only RepliSeq)
    unite("mb_domain", contains(chromatin_features$name), sep = "_") %>% 
    # binarize weighted average DNA repair value by being lower or larger than the across-genome median
    rowwise %>% 
    mutate_at(vars(contains(match = dnarep_marks$name)),
              function(x){var_name = rlang::as_label(substitute(x))
              ifelse(!is.na(suppressWarnings(as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score))),
                     # it's numeric score, binarize as low/high
                     yes = ifelse(x <= as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score),
                                  yes = "low",
                                  no = "high"),
                     # it's factor, leave as is
                     no = x)}) %>% 
    # dna rep mark levels as ordered factors
    mutate_at(vars(contains(match = dnarep_marks$name)),
              ~if(unique(.)[1] %in% c('AID_target', 'bgGenome')){
                factor(., ordered = T, levels = c('AID_target', 'bgGenome')) # higher mut rates --> baseline
              }else{
                factor(., ordered = T, levels = c('low', 'high'))}) # baseline --> lower mut rates
  
  list_failed_samples[[sample]] = merged %>%
    select(-mut_id) %>% 
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
  rm(merged)
  gc()
}

bound_failed_samples = list_failed_samples %>% 
  map2_df(., names(.), ~mutate(.x, sample = .y)) %>% 
  bind_rows %>% 
  mutate(sample_type = "10 failed samples") %>% 
  relocate(sample, sample_type)


######################################################
## same for 20 first successful samples

first20_successful_samples = read_tsv("../successful", col_names = F)$X1

list_first20_successful_samples = list()
gc()

for(sample in first20_successful_samples){
  
  existing_file = c()
  for(file_exists in paste0(path_somatic_variation, sample, ".csv")){
    if(file.exists(file_exists)){
      existing_file = c(existing_file, file_exists)
    }
  }
  if(length(existing_file) != 1){
    stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
  }
  
  dfright = read_csv(existing_file) %>%
    select(chr, start, end, tri) %>% 
    rename("chrom" = "chr") %>%
    mutate(mut_id = paste0("mut_", row_number()))
  
  merged = dfright %>%
    bed_intersect(dfleft, suffix = c("_dfright", "_dfleft")) %>%
    select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
    # remove the "dleft" and "dright" parts of the column names
    rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
    # combine chromatin features (although there should typically be only RepliSeq)
    unite("mb_domain", contains(chromatin_features$name), sep = "_") %>% 
    # binarize weighted average DNA repair value by being lower or larger than the across-genome median
    rowwise %>% 
    mutate_at(vars(contains(match = dnarep_marks$name)),
              function(x){var_name = rlang::as_label(substitute(x))
              ifelse(!is.na(suppressWarnings(as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score))),
                     # it's numeric score, binarize as low/high
                     yes = ifelse(x <= as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score),
                                  yes = "low",
                                  no = "high"),
                     # it's factor, leave as is
                     no = x)}) %>% 
    # dna rep mark levels as ordered factors
    mutate_at(vars(contains(match = dnarep_marks$name)),
              ~if(unique(.)[1] %in% c('AID_target', 'bgGenome')){
                factor(., ordered = T, levels = c('AID_target', 'bgGenome')) # higher mut rates --> baseline
              }else{
                factor(., ordered = T, levels = c('low', 'high'))}) # baseline --> lower mut rates
  
  list_first20_successful_samples[[sample]] = merged %>%
    select(-mut_id) %>% 
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
  rm(merged)
  gc()
}

bound_first20_successful_samples = list_first20_successful_samples %>% 
  map2_df(., names(.), ~mutate(.x, sample = .y)) %>% 
  bind_rows %>% 
  mutate(sample_type = "first 20 successful samples") %>% 
  relocate(sample, sample_type)



####################################################################
### plot

bound_samples = bind_rows(bound_first20_successful_samples, bound_failed_samples) %>% 
  mutate(sample = gsub("TCGA-FF-", "TCGA_FF_", sample),
         sample = gsub("-.*", "", sample))
gc()

snvs_per_bin_plot = ggplot(bound_samples,
       aes(x = mutcount,
           fill = sample)) +
  geom_bar(width = 0.6,
           position = position_dodge()) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  facet_grid(facets = ~sample_type,
             space = "free",
             scales = "free") +
  theme_classic() +
  xlab("#SNVs at a genomic bin (sqrt scale)") +
  ylab("#Bins (10 binary features × 6 RepliSeq × 96 SBS ; sqrt scale)") +
  theme(text = element_text(size = 12),
        legend.position = "none")
ggsave("snvs_per_bin_plot.jpg",
       plot = snvs_per_bin_plot,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "white")
