library(tidyverse)
library(dtplyr)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19") # all coordinates are in hg19
library(rlang)
library(spgs)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")


# from command parameters
args = commandArgs(trailingOnly=TRUE)

# load utils.R (functions) -- only used here if trinucleotide matching
if(interactive()){
  source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")
} else {
  source(args[1])
}

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[2]) %>%
  read_csv(comment = "#")

# load raw-score feature map from previous process
map_features = ifelse(interactive(),
                      yes = "../work/67/2ffa9d1db5e0ed41abf08e1190cb1e/map_features_single_chr",
                      no = args[3]) %>% 
  read_tsv()

chromosome = unique(map_features$seqnames)


# load collected median_scores from 1st process
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
                       no = lapply(list(args[-(1:3)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)


#### binarize scores, and get frequency of each trinucleotide per sequence (for offset process)

feature_names = map_features %>% 
  select(-c("seqnames", "start", "end", "width", "strand")) %>% 
  names
features_with_numeric_score = map_features %>% 
  select(-c("seqnames", "start", "end", "width", "strand")) %>% 
  select_if(is.numeric) %>% 
  names
features_with_character_levels = map_features %>% 
  select(-c("seqnames", "start", "end", "width", "strand")) %>% 
  select_if(is.character) %>% 
  names

## binarize weighted average DNA repair value by being lower or larger than the across-genome median

map_features_binarized_temp = map_features %>%
  lazy_dt %>% 
  #### WARNING first do the average score at duplicated (start end) ranges, this is due to the (in some features) hg38-->hg19 lift dividing some ranges into 2 alternative ranges with the same score
  group_by_at(vars('seqnames', 'start', 'end', features_with_character_levels)) %>% 
  summarise_at(features_with_numeric_score,
               ~mean(.)) %>% 
  ungroup %>% 
  as_tibble %>% 
  rowwise %>% 
  lazy_dt %>% 
  mutate_at(vars(features_with_numeric_score[!features_with_numeric_score %in% chromatin_features$name]),
            function(x){var_name = rlang::as_label(substitute(x))
            ifelse(x <= filter(median_scores, dnarep_mark == var_name) %>% pull(median_score),
                   "low",
                   "high")}) %>% 
  as_tibble %>% 
  relocate(features_with_character_levels, .after = last_col()) %>% 
  unite("metadata", !matches("seqnames|start|end|width|strand")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

## collapse contiguous ranges if they have same metadata levels
map_features_binarized_temp = unlist(reduce(split(map_features_binarized_temp, ~metadata)))
mcols(map_features_binarized_temp) = names(map_features_binarized_temp)
map_features_binarized = map_features_binarized_temp %>% 
  as_tibble %>% 
  arrange(start) %>% 
  mutate(X = gsub("AID_", "AID ", X)) %>% 
  separate(X, into = feature_names, sep = "_") %>% 
  mutate(AID_regions = gsub("AID ", "AID_", AID_regions)) %>% 
  # add 1 bp downstream and upstream regions of width 1 or 2, so trinucs can be fetched
  mutate(start = ifelse(width <= 2,
                        start-1,
                        start),
         end = ifelse(width <= 2,
                      end+1,
                      end),
         width = end-start+1)
gc()


# get sequence for each range
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                   names = makeGRangesFromDataFrame(map_features_binarized, 
                                                    keep.extra.columns = T))


# get frequency of each trinucleotide per sequence, moving 1 nucleotide downstream each time
trinuc32_freq = trinucleotideFrequency(sequences) %>%
  as_tibble %>%
  rownames_to_column("id") %>% 
  lazy_dt %>%
  pivot_longer(cols = -id,
               names_to = 'trinuc32',
               values_to = "freq") %>%
  group_by(id) %>%
  # trinucs that do not have a C or T in center, convert to reverse complement
  mutate(trinuc32 = ifelse(substr(trinuc32, start = 2, stop = 2) %in% c('A', 'G'),
                           spgs::reverseComplement(trinuc32, case="upper"),
                           trinuc32)) %>% 
  # sum up frequencies of each N(C|T)N & reverse complement pair, within id
  group_by(id, trinuc32) %>%
  summarise(freq = sum(freq)) %>%
  ungroup %>% 
  # back to original format (each row ('id') maps to the same row in map_features)
  pivot_wider(names_from = 'trinuc32',
              values_from = 'freq') %>% 
  as_tibble %>% 
  arrange(as.numeric(id)) %>% 
  select(-id)
gc()

# bind trinuc32 freqs to map_features_binarized
map_features_binarized_trinuc32_freq = map_features_binarized %>% 
  bind_cols(trinuc32_freq) %>% 
  group_by_at(vars(feature_names)) %>% 
  summarise_at(vars(matches("^[A,C,T,G][C,T][A,C,T,G]$")),
               ~sum(.)) %>% 
  mutate(chrom = chromosome) %>% 
  relocate(chrom)
gc()

write_tsv(map_features_binarized_trinuc32_freq, paste0("map_features_binarized_", chromosome, ".tsv"))
