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


### load sample; merge with dna repair, chromatin landscape, and offset; parse; regression

# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample = ifelse(interactive(),
                yes = "MSM0.128", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_",
                                no = args[2]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

metadata_sample = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv",
                         no = args[3]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info"))) %>% 
  filter(sample_id == sample)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[4]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[5]) %>%
  read_csv(comment = "#")


# load offset from 3rd process
offset = ifelse(interactive(),
                yes = "offset.tsv",
                no = args[6]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"


## load map_features (all chromosomes) from 2nd process
dfleft = ifelse(interactive(),
                yes = "map_features_chr21.tsv", #"../res/map_features.tsv"
                no = paste0(args[7], "/res/map_features.tsv")) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()


## DONT filter out low mappability regions
# low_mappability_regions = ifelse(interactive(),
#                                  yes = "/g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed",
#                                  no = args[8]) %>% 
#   import.bed() %>% data.frame %>%
#   mutate(seqnames = gsub("^", "chr", seqnames)) %>% 
#   rename("chrom" = "seqnames")


# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/26/08b67aa330177901dc72e70df9cc6b/median_score_OGG1_GOx30_chipseq.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/02/bd5527ea98ae2fbcb50d1c7947a0fb/median_score_OGG1_GOx60_chipseq.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/fb/0c57088546769b3d33c5ce901bdda8/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/f2/acd0d55a2fb0d7a07260a120f0d01e/median_score_UV_XRseq_NHF1_CPD_1h.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/d2/65d1e94a5015f7106d266899bae572/median_score_XRCC4.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/58/2192db8e29c53114153463b867ce68/median_score_SETD2_control.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/73/4d3c3f1693923f3808026f1194d58f/median_score_MSH6_control.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/35/ba92f745984eec59ec4ca62cd3d2df/median_score_TP53_dauno_MOLM13.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/db/db4afaf45c9ec6f18e486dc5dfcc87/median_score_TP53_dauno_K562.tsv",
                                           "../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/fb/8caddf7603a8d872bd5df497010bea/median_score_AID_regions.tsv")), 
                                    read_tsv),
                       no = lapply(list(args[-(1:8)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)


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

# load it
dfright = read_csv(existing_file) %>%
  select(chr, start, end) %>% 
  rename("chrom" = "chr") %>%
  mutate(mut_id = paste0("mut_", row_number()))
gc()


### map chromatin features
merged = dfright %>%
  # DONT remove low mappability regions from SNVs table
  #bed_subtract(low_mappability_regions) %>% 
  bed_intersect(dfleft, suffix = c("_dfright", "_dfleft")) %>%
  select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"), contains("overlap"))) %>%
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
gc()

if(nrow(merged) == 0){
  stop("ERROR - Empty 'merged' table: probably 'low_mappability_regions' has removed all SNVs for this sample! Exiting...\n")
}


## add mutation counts per mb_domain
sommut_dnarep_chromatin = merged %>%
  select(-mut_id) %>% 
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

write_tsv(results_sample, "results_sample.tsv")
