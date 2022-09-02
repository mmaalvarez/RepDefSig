library(tidyverse)
library(data.table)
# BiocManager::install('XVector', 'GenomicRanges', 'rtracklayer')
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
library(MASS)
library(lme4)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")


### load sample; merge with dna repair, chromatin landscape, and offset; parse; regression


# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample = ifelse(interactive(),
                yes = "CPCT02010350T",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_",
                                no = args[2])

metadata_sample = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv",
                         no = args[3]) %>%
  read_tsv %>% filter(sample_id == sample)

dnarep_marks = ifelse(interactive(),
                      yes = "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/1_parser_and_regressions/model_0_test/input_lists/dnarep_marks.tsv",
                      no = args[4]) %>%
  read_tsv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/1_parser_and_regressions/model_0_test/input_lists/chromatin_features.tsv",
                            no = args[5]) %>%
  read_tsv(comment = "#")


# load results from previous process

ln_bpsum_chromatin_env_table = read_tsv("ln_bpsum_chromatin_env_table.tsv")

map_features = as_tibble(fread("map_features.tsv"))
gc()


## load sample
somatic_mutations_granges = read_csv(paste0(path_somatic_variation, sample, ".csv")) %>%
  select(-c(sample_id, ref, alt, alteration, context)) %>%
  # add +-50Kb buffer around each SNV (i.e. 1/10Mb windows), to ensure most of them overlap with at least 1 reptime genomic range
  mutate(start = start - 50000,
         end = end + 50000) %>%
  makeGRangesFromDataFrame(keep.extra.columns=T)


## map chromatin features
dfleft = data.frame(map_features) %>% 
  rename("chrom" = "seqnames")

dfright = data.frame(somatic_mutations_granges) %>% 
  rename("chrom" = "seqnames") %>%
  mutate(mut_id = paste0("mut_", row_number()))

merged = bed_intersect(dfleft, dfright, suffix = c("_dfleft",
                                                   "_dfright")) %>%
  select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
  ## weighted averages of chromatin feature bins and dna repair abundances
  mutate(`.overlap` = `.overlap` + 1) %>%
  group_by(mut_id_dfright) %>%
  summarise(total_overlaps_length = sum(`.overlap`), across()) %>%
  mutate_at(vars(!contains("mut_id_dfright") & !contains("total_overlaps_length") & !contains("tri_dfright") & !contains(".overlap")),
            ~. * `.overlap` / total_overlaps_length) %>%
  group_by(mut_id_dfright, tri_dfright) %>%
  summarise_at(vars(!contains("mut_id_dfright") & !contains("total_overlaps_length") & !contains("tri_dfright") & !contains(".overlap")),
               ~sum(.)) %>%
  # round the reptime bins to 0 decimal
  mutate_at(vars(contains(chromatin_features$name)),
            ~round(., 0))
colnames(merged) = gsub("_dfleft", "", colnames(merged))
colnames(merged) = gsub("_dfright", "", colnames(merged))

# combine chromatin features (although there should typically be only RepliSeq)
merged = merged %>%
  unite("chromatin_env", contains(chromatin_features$name), sep = "_")

gc()


## add mutation counts per trinuc×chromenv
sommut_tricount_dnarep_chromatin = table(merged$chromatin_env,
                                         merged$tri) %>%
  as.data.frame() %>%
  rename("chromatin_env" = "Var1",
         "tri" = "Var2",
         "mutcount" = "Freq") %>%
  # dont care about trinuc×chromenv with no reported SNVs for this sample
  filter(mutcount != 0) %>%
  # append dna repair abundances
  merge(merged, all = T) %>%
  # average dna repair abundances per chromatin_env×trinuc pair that has >=2 mutations
  group_by(tri, chromatin_env, mutcount) %>%
  summarise_at(vars(!contains("tri") & !contains("chromatin_env") & !contains("mutcount") & !contains("mut_id")),
               ~mean(.)) %>% 
  ungroup
# parse
sommut_tricount_dnarep_chromatin$chromatin_env = factor(sommut_tricount_dnarep_chromatin$chromatin_env, ordered = T)
trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")
sommut_tricount_dnarep_chromatin$tri = factor(sommut_tricount_dnarep_chromatin$tri, ordered = T, levels = trinuc_96)
sommut_tricount_dnarep_chromatin = sommut_tricount_dnarep_chromatin %>%
  ## add offset
  merge(ln_bpsum_chromatin_env_table) %>%
  relocate(mutcount) %>%
  relocate(chromatin_env, .before = "ln_bpsum_chromatin_env") %>%
  relocate(tri, .after = "chromatin_env") %>%
  arrange(tri, chromatin_env)

rm(merged) ; gc()


### regression

formula = paste0("mutcount ~ ",
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "(1 | chromatin_env) + ",
                 "(1 | tri) + ",
                 "offset(ln_bpsum_chromatin_env)")

# first try a generalized linear mixed-effects model for the negative binomial family
y = tryCatch(glmer.nb(formula = formula, 
                      data = sommut_tricount_dnarep_chromatin),
             # if there is a failure to converge to theta, run Poisson
             warning = function(w) glmer(formula = formula, 
                                         data = sommut_tricount_dnarep_chromatin, 
                                         family = poisson),
             error = function(e) glmer(formula = formula, 
                                       data = sommut_tricount_dnarep_chromatin, 
                                       family = poisson)) %>%
  # parse output
  broom.mixed::tidy() %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  select(term, estimate, p.value) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>%
  mutate(sample_id = sample)


## append features' coefficients and pvalues to metadata_sample
results_sample = full_join(metadata_sample, y) %>%
  relocate(source) %>% relocate(sample_id_2) %>% relocate(sample_id)

write_tsv(results_sample, "results_sample.tsv")
gc()
