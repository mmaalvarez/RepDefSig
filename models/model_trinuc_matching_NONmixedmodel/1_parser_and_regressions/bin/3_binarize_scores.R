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


# load utils.R (functions) -- only used here if trinucleotide matching
if(interactive()){
  source("../../../../bin/utils.R")
} else {
  source(args[1])
}

# from command parameters
args = commandArgs(trailingOnly=TRUE)

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[1]) %>%
  read_csv(comment = "#")


## load raw-score feature map from previous process
if(interactive()){ ## if interactive, load all chromosomes obtained in the non-trinuc-matching model, as they can be reused in this trinuc-matching model
  
  # retrieve them
  map_features_other_model = c(
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr1.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr2.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr3.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr4.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr5.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr6.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr7.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr8.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr9.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr10.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr11.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr12.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr13.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr14.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr15.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr16.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr17.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr18.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr19.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr20.tsv")),
    Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr21.tsv"))#,
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chr22.tsv")),
    # Sys.glob(paste0("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/map_features_chrX.tsv"))
    )
  
  # load them
  map_features_other_model_dataFiles = lapply(map_features_other_model,
                                              read_tsv) ; gc()
  # name them by their chromosome
  for(dataFile in seq(1, length(map_features_other_model_dataFiles), 1)){
    names(map_features_other_model_dataFiles)[[dataFile]] = unique(map_features_other_model_dataFiles[[dataFile]]$seqnames) ; gc()
    } ; gc()
  
  # choose one chromosome
  chromosome = "chr21"
  map_features = map_features_other_model_dataFiles[[chromosome]]
  
  } else {
    # if interactive, load just one chromosome at a time (in parallel with nextflow), from this model
    map_features = read_tsv(args[2])
    
    chromosome = unique(map_features$seqnames)
  }


# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c(Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_OGG1_GOx30_chipseq.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_OGG1_GOx60_chipseq.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_UV_XRseq_NHF1_CPD_1h.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_XRCC4.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_SETD2_control.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_MSH6_control.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_TP53_dauno_K562.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_TP53_dauno_MOLM13.tsv")[1],
                                           Sys.glob("../../../model_NONmatching_mixedmodel/1_parser_and_regressions/work/[[:alnum:]][[:alnum:]]/*/median_score_AID_regions.tsv")[1])), 
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

map_features_binarized = map_features %>%
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
rm(map_features) ; gc()

## collapse contiguous ranges if they have same metadata levels
map_features_binarized = unlist(reduce(split(map_features_binarized, ~metadata)))
mcols(map_features_binarized) = names(map_features_binarized)
map_features_binarized = map_features_binarized %>% 
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
rm(sequences) ; gc()

# bind trinuc32 freqs to map_features_binarized
map_features_binarized_trinuc32_freq = map_features_binarized %>% 
  bind_cols(trinuc32_freq) %>% 
  group_by_at(vars(feature_names)) %>% 
  summarise_at(vars(matches("^[A,C,T,G][C,T][A,C,T,G]$")),
               ~sum(.)) %>% 
  mutate(chrom = chromosome) %>% 
  relocate(chrom) %>% 
  # prepare for matching
  unite(col = "bin", !matches("^[A,C,T,G][C,T][A,C,T,G]$")) %>% 
  column_to_rownames("bin") %>% 
  # remove "all-0-mut" rows (bins)
  filter(rowSums(.) >= 1)
rm(trinuc32_freq) ; rm(map_features_binarized) ; gc()

### trinuc matching
map_features_binarized_trinuc32_freq_matched = trinuc_matching(map_features_binarized_trinuc32_freq,
                                                               stoppingCriterion = 0.001,
                                                               maxIter = 20000*length(map_features_binarized_trinuc32_freq),
                                                               n_finish_tokens = 1000,
                                                               max_fraction_removed_muts = 0.25) %>% 
  mutate(bin = gsub("AID_", "AID", bin)) %>% 
  separate(bin, into = map_features_binarized %>% 
                         select(-c(start, end, width, strand)) %>% 
                         names) %>% 
  mutate_all(~gsub("AIDtarget", "AID_target", .))
rm(map_features_binarized_trinuc32_freq) ; gc()

write_tsv(map_features_binarized_trinuc32_freq_matched, paste0("map_features_binarized_", chromosome, ".tsv"))
