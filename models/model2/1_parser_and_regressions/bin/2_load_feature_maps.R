library(tidyverse)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
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


### load chromatin feature maps and dna repair marks hg19 and keep as granges

# from command parameters
args = commandArgs(trailingOnly=TRUE)

dnarep_marks = ifelse(interactive(),
                      yes = "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model2/1_parser_and_regressions/work/ad/22931b2720f1cfd778be28079a1f02/dnarep_marks", #"../input_lists/dnarep_marks.csv",
                      no = args[1]) %>%
  read_csv(comment = "#") 

chromatin_features = ifelse(interactive(),
                            yes = "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model2/1_parser_and_regressions/work/ad/22931b2720f1cfd778be28079a1f02/chromatin_features", #"../input_lists/chromatin_features.csv",
                            no = args[2]) %>%
  read_csv(comment = "#")

chromosome = ifelse(interactive(),
                    yes = "chr9",
                    no = paste0("chr", args[3]))

low_mappability_regions = read_tsv(args[4], col_names = F)

# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c("../work/ad/22931b2720f1cfd778be28079a1f02/median_score_OGG1_GOx30_chipseq.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_MSH6_control.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_OGG1_GOx60_chipseq.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_SETD2_control.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_UV_XRseq_NHF1_CPD_1h.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv",
                                           "../work/ad/22931b2720f1cfd778be28079a1f02/median_score_XRCC4.tsv",
                                           "../work/XXX/XXX/median_score_TP53_dauno_K562.tsv",
                                           "../work/XXX/XXX/median_score_TP53_dauno_MOLM13.tsv")), 
                                    read_tsv),
                       no = lapply(list(args[-(1:4)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)


dnarep_marks_list_files = list()
for (feature in dnarep_marks$name){ #[1]){

  path_file = filter(dnarep_marks, name == feature)$path

  dnarep_marks_list_files[[feature]] = tryCatch(import.bw(path_file),
                                                error = function(e) tryCatch(import.bedGraph(path_file),
                                                                             error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file),
                                                                                                                                   keep.extra.columns = T),
                                                                                                          error = function(e) import.bed(path_file))))
  # extract chromosome (args[3]), otherwise it gets too long and the bed_intersect crashes
  dnarep_marks_list_files[[feature]] = dnarep_marks_list_files[[feature]][seqnames(dnarep_marks_list_files[[feature]]) == chromosome]

  # add feature name as the metadata (score) colname
  colnames(elementMetadata(dnarep_marks_list_files[[feature]])) = feature
  gc()
}

chromatin_features_list_files = list()
for (feature in chromatin_features$name){

  path_file = filter(chromatin_features, name == feature)$path

  chromatin_features_list_files[[feature]] = tryCatch(import.bw(path_file),
                                                      error = function(e) tryCatch(import.bedGraph(path_file),
                                                                                   error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file),
                                                                                                                                         keep.extra.columns = T),
                                                                                                                error = function(e) import.bed(path_file))))
  chromatin_features_list_files[[feature]] = chromatin_features_list_files[[feature]][seqnames(chromatin_features_list_files[[feature]]) == chromosome]
  gc()
}


## merge Reptime and dna repair coordinates

dfleft = data.frame(chromatin_features_list_files[[1]]) %>% rename("chrom" = "seqnames")

list_feature_names = list()
list_feature_names[[1]] = colnames(elementMetadata(chromatin_features_list_files[[1]]))
n_chromatin_features = length(unlist(list_feature_names))

for (feature_i in seq(1, length(dnarep_marks_list_files))){

  dfright = data.frame(dnarep_marks_list_files[[feature_i]]) %>% 
    rename("chrom" = "seqnames")
  
  list_feature_names[[n_chromatin_features + feature_i]] = colnames(elementMetadata(dnarep_marks_list_files[[feature_i]]))

  merged_temp = bed_intersect(dfleft, dfright, suffix = c("_dfleft",
                                                          "_dfright")) %>%
    mutate(start = ifelse(start_dfleft >= start_dfright,
                          start_dfleft,
                          start_dfright),
           end = ifelse(end_dfleft <= end_dfright,
                        end_dfleft,
                        end_dfright),
           strand = "*",
           width = end - start + 1) %>%
    select(chrom, start, end, width, strand, contains(unlist(list_feature_names)))
  colnames(merged_temp) = gsub("_dfleft", "", colnames(merged_temp))
  colnames(merged_temp) = gsub("_dfright", "", colnames(merged_temp))

  dfleft = merged_temp
  rm(merged_temp) ; rm(dfright)
  gc()
}

rm(dnarep_marks_list_files) ; gc()

map_features = dfleft %>%
  rename("seqnames" = "chrom") %>%
  arrange(start)

rm(dfleft) ; gc()

# this will go to last process
write_tsv(map_features, paste0("map_features_", chromosome, ".tsv"))
# map_features = read_tsv("map_features_chr21.tsv")


#### binarize scores, and get frequency of each trinucleotide per sequence (for 3rd process -- offset)

## binarize weighted average DNA repair value by being lower or larger than the across-genome median
map_features_binarized_temp = map_features %>%
  lazy_dt %>% 
  #### WARNING first do the average score at duplicated (start end) ranges, this is due to the hg38-->hg19 lift dividing some ranges into 2 alternative ranges with the same score
  group_by(seqnames, start, end) %>% 
  summarise_at(vars(!matches("seqnames|start|end|width|strand")),
               ~mean(.)) %>% 
  # there are some duplicated only at the start
  group_by(seqnames, start) %>% 
  summarise_at(vars(!matches("seqnames|start|width|strand")),
               ~mean(.)) %>% 
  # there are some duplicated only at the start
  group_by(seqnames, end) %>% 
  summarise_at(vars(!matches("seqnames|end|width|strand")),
               ~mean(.)) %>% 
  ungroup %>% 
  as_tibble %>% 
  rowwise %>% 
  lazy_dt %>% 
  mutate_at(vars(dnarep_marks$name),
            function(x){var_name = rlang::as_label(substitute(x))
            ifelse(x <= filter(median_scores, dnarep_mark == var_name) %>% pull(median_score),
                   "low",
                   "high")}) %>% 
  as_tibble %>% 
  unite("metadata", !matches("seqnames|start|end|width|strand")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
## collapse contiguous ranges if they have same metadata levels
map_features_binarized_temp = unlist(reduce(split(map_features_binarized_temp, ~metadata)))
mcols(map_features_binarized_temp) = names(map_features_binarized_temp)
map_features_binarized = map_features_binarized_temp %>% 
  as_tibble %>% 
  arrange(start) %>% 
  separate(X, into = c(chromatin_features$name, dnarep_marks$name), sep = "_") %>% 
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


### TRINUC MATCHING


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
  group_by_at(vars(c(chromatin_features$name, dnarep_marks$name))) %>% 
  summarise_at(vars(matches("^[A,C,T,G][C,T][A,C,T,G]$")),
               ~sum(.)) %>% 
  mutate(chrom = chromosome) %>% 
  relocate(chrom)
gc()

write_tsv(map_features_binarized_trinuc32_freq, paste0("map_features_binarized_", chromosome, ".tsv"))
