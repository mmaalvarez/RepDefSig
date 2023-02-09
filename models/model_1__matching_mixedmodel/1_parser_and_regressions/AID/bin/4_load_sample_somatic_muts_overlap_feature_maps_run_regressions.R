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


### load sample; merge with AID regions, chromatin landscape, and offset; parse; regression

# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample = ifelse(interactive(),
                yes = "0206c56a-d7f2-5800-8444-3df54debd7a0",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_",
                                no = args[2]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)


AID_regions = ifelse(interactive(),
                      yes = "../input_lists/AID_regions.csv",
                      no = args[3]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[4]) %>%
  read_csv(comment = "#")


# load offset from previous process
offset = ifelse(interactive(),
                yes = "./offset.tsv",
                no = args[5]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"


## load map_features (all chromosomes) from 1st process
dfleft = ifelse(interactive(),
                yes = "map_features_chr21.tsv",
                no = paste0(args[6], "/res/map_features.tsv")) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()


low_mappability_regions = read_tsv(args[7], col_names = F)


## load sample (somatic mutations)

# find path+sample name that exists
existing_file = c()
for(file_exists in paste0(path_somatic_variation, sample, ".csv")){
  if(file.exists(file_exists)){
    existing_file = c(existing_file, file_exists)
  }
}

# raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
if(length(existing_file) != 1){
  stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
}

# load
dfright = read_csv(existing_file) %>%
  select(chr, start, end, tri) %>% 
  rename("chrom" = "chr") %>%
  mutate(mut_id = paste0("mut_", row_number()))
gc()


### map chromatin features
merged = bed_intersect(dfleft, dfright, suffix = c("_dfleft",
                                                   "_dfright")) %>%
  select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
  # remove the "dleft" and "dright" parts of the column names
  rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
  # combine chromatin features (although there should typically be only RepliSeq)
  unite("mb_domain", contains(chromatin_features$name), sep = "_")

gc()


## add mutation counts per trinuc×mb_domain
trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")
sommut_tricount_AIDregions_chromatin = merged %>%
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
  # dna rep mark levels as ordered factors
  mutate_at(vars(contains(match = AID_regions$name)),
            ~ factor(., ordered = T, levels = c('bgGenome', 'AID_target'))) %>% 
  # mb_domain and tri as ordered and unordered factors, respectively
  mutate(mb_domain = factor(mb_domain, ordered = T),
         tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
  arrange(tri, mb_domain) %>% 
  as_tibble

rm(merged) ; gc()


### regression

formula = paste0("mutcount ~ ",
                 paste(AID_regions$name, collapse = " + "), " + ",
                 "(1 | mb_domain) + ",
                 "(1 | tri) + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear mixed-effects model for the negative binomial family
y = suppressWarnings(glmer.nb(formula = formula, 
                     data = sommut_tricount_AIDregions_chromatin))

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


write_tsv(y_tidy, "results_sample.tsv")
