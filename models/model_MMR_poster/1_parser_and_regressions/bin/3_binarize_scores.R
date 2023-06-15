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


chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[1]) %>%
  read_csv(comment = "#")


## load raw-score feature map from previous process
map_features = ifelse(interactive(),
                      yes = Sys.glob("map_features_chr21.tsv"),
                      no = args[2]) %>% 
  read_tsv()

chromosome = unique(map_features$seqnames)


# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c(Sys.glob("median_score_DHS.tsv")[1],
                                           Sys.glob("median_score_exons.tsv")[1],
                                           Sys.glob("median_score_H3K36me3.tsv")[1],
                                           Sys.glob("median_score_RepliSeq.tsv")[1],
                                           Sys.glob("median_score_RnaSeq.tsv")[1])), 
                                    read_tsv),
                       no = lapply(list(args[-(1:2)]), read_tsv)) %>%
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
  group_by_at(vars('seqnames', 'start', 'end', all_of(features_with_character_levels))) %>% 
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
  relocate(all_of(features_with_character_levels), .after = last_col()) %>% 
  unite("metadata", !matches("seqnames|start|end|width|strand")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
rm(map_features) ; gc()

## collapse contiguous ranges if they have same metadata levels
map_features_binarized_temp = unlist(reduce(split(map_features_binarized_temp, ~metadata)))
mcols(map_features_binarized_temp) = names(map_features_binarized_temp)
map_features_binarized = map_features_binarized_temp %>% 
  as_tibble %>% 
  arrange(start) %>% 
  mutate(#X = gsub("CTCF_cohesin_peak", "CTCFcohesinpeak ", X),
         #X = gsub("hairpin_TpCpH", "hairpinTpCpH ", X),
         X = gsub("exon__non_CEG_TSG_OG", "exon", X)) %>%
  separate(X, into = feature_names, sep = "_") %>% 
  mutate(#CTCF_cohesin = gsub("CTCFcohesinpeak ", "CTCF_cohesin_peak", CTCF_cohesin),
         #A3A_TpCpH_hairpins = gsub("hairpinTpCpH ", "hairpin_TpCpH", A3A_TpCpH_hairpins),
         exons = gsub("exon", "exon__non_CEG_TSG_OG", exons)) %>% 
  # add 1 bp downstream and upstream regions of width 1 or 2, so trinucs can be fetched
  mutate(start = ifelse(width <= 2,
                        start-1,
                        start),
         end = ifelse(width <= 2,
                      end+1,
                      end),
         width = end-start+1)
rm(map_features_binarized_temp) ; gc()


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
rm(sequences) ; gc()

# bind trinuc32 freqs to map_features_binarized
map_features_binarized_trinuc32_freq = map_features_binarized %>% 
  bind_cols(trinuc32_freq) %>% 
  group_by_at(vars(all_of(feature_names))) %>% 
  summarise_at(vars(matches("^[A,C,T,G][C,T][A,C,T,G]$")),
               ~sum(.)) %>% 
  mutate(chrom = chromosome) %>% 
  relocate(chrom)

rm(trinuc32_freq) ; rm(map_features_binarized) ; gc()

write_tsv(map_features_binarized_trinuc32_freq, paste0("map_features_binarized_", chromosome, ".tsv"))
