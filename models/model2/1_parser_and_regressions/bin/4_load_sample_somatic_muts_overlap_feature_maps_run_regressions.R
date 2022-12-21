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
                yes = "MSM0.1", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_",
                                no = args[2]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

metadata_sample = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv",
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
                yes = "map_features_chr21.tsv",
                no = paste0(args[7], "/res/map_features.tsv")) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()


# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply("median_score_OGG1_GOx30_chipseq.tsv", read_tsv),
                       no = lapply(args[-(1:7)], read_tsv)) %>%
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

# load
dfright = read_csv(existing_file) %>%
  select(chr, start, end, tri) %>% 
  # add +-50bp buffer around each SNV (i.e. 0.1kb windows), to ensure most of them overlap with reptime and DNA mark genomic ranges (I still lose quite a lot)
  mutate(start = start - 50,
         end = end + 50) %>%
  rename("chrom" = "chr") %>%
  mutate(mut_id = paste0("mut_", row_number()))
gc()


### map chromatin features
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
                               "high")})
gc()


## add mutation counts per trinuc×mb_domain
trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")
sommut_tricount_dnarep_chromatin = merged %>%
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
  # dna rep mark levels as ordered factors
  mutate_at(vars(contains(match = dnarep_marks$name)),
            ~ factor(., ordered = T, levels = c('low', 'high'))) %>% 
  # mb_domain and tri as ordered and unordered factors, respectively
  mutate(mb_domain = factor(mb_domain, ordered = T),
         tri = factor(tri, ordered = F, levels = trinuc_96)) %>% 
  arrange(tri, mb_domain) %>% 
  as_tibble

rm(merged) ; gc()


### regression

formula = paste0("mutcount ~ ",
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "(1 | mb_domain) + ",
                 "(1 | tri) + ",
                 "offset(log_freq_trinuc32)")

# first try a generalized linear mixed-effects model for the negative binomial family
y = tryCatch(glmer.nb(formula = formula, 
                      data = sommut_tricount_dnarep_chromatin),
             # if there is a failure to converge to theta, run Poisson
             warning = function(w) glmer(formula = formula, 
                                         data = sommut_tricount_dnarep_chromatin, 
                                         family = poisson),
             error = function(e) glmer(formula = formula, 
                                       data = sommut_tricount_dnarep_chromatin, 
                                       family = poisson))
# track whether NB or Poisson
nb_or_pois = family(y)[[1]]

# parse output (calc 95%CI same as with stats::confint())
y = broom.mixed::tidy(y,
                      conf.int = T, conf.method = "profile", exponentiate = F, effects = "fixed") %>% 
  filter(effect == "fixed" & term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high, p.value)) %>%
  mutate(sample_id = sample,
         glm = ifelse(str_detect(nb_or_pois, "Negative Binomial"),
                      "Negative Binomial",
                      ifelse(str_detect(nb_or_pois, "[P,p]oisson"),
                             "Poisson",
                             stop(paste0("ERROR! GLM family used was not 'Negative Binomial' nor '[P,p]oisson'; instead, it was '", nb_or_pois, "'")))))
gc()

## append features' coefficients and pvalues to metadata_sample
results_sample = full_join(metadata_sample, y) %>%
  relocate(sample_id)

write_tsv(results_sample, "results_sample.tsv")
