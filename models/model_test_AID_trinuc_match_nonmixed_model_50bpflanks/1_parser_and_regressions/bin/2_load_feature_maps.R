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


### load chromatin feature maps and AID regions marks hg19 and keep as granges

# from command parameters
args = commandArgs(trailingOnly=TRUE)

AID_regions = ifelse(interactive(),
                      yes = "../input_lists/AID_regions.csv",
                      no = args[1]) %>%
  read_csv(comment = "#") 

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[2]) %>%
  read_csv(comment = "#")

chromosome = ifelse(interactive(),
                    yes = "chr21",
                    no = paste0("chr", args[3]))

low_mappability_regions = read_tsv(args[4], col_names = F)


AID_regions_list_files = list()
for (feature in AID_regions$name){

  path_file = filter(AID_regions, name == feature)$path

  AID_regions_list_files[[feature]] = tryCatch(import.bw(path_file),
                                                error = function(e) tryCatch(import.bedGraph(path_file),
                                                                             error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file),
                                                                                                                                   keep.extra.columns = T),
                                                                                                          error = function(e) import.bed(path_file))))
  # extract chromosome (args[3]), otherwise it gets too long and the bed_intersect crashes
  AID_regions_list_files[[feature]] = AID_regions_list_files[[feature]][seqnames(AID_regions_list_files[[feature]]) == chromosome]

  # add feature name as the metadata (score) colname
  colnames(elementMetadata(AID_regions_list_files[[feature]])) = feature
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


## merge Reptime and AID regions coordinates

dfleft = data.frame(chromatin_features_list_files[[1]]) %>% rename("chrom" = "seqnames")

list_feature_names = list()
list_feature_names[[1]] = colnames(elementMetadata(chromatin_features_list_files[[1]]))
n_chromatin_features = length(unlist(list_feature_names))

for (feature_i in seq(1, length(AID_regions_list_files))){

  dfright = data.frame(AID_regions_list_files[[feature_i]]) %>% 
    rename("chrom" = "seqnames")
  
  list_feature_names[[n_chromatin_features + feature_i]] = colnames(elementMetadata(AID_regions_list_files[[feature_i]]))

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

rm(AID_regions_list_files) ; gc()

map_features = dfleft %>%
  rename("seqnames" = "chrom") %>%
  arrange(start)

rm(dfleft) ; gc()

# this will go to last process
write_tsv(map_features, paste0("map_features_", chromosome, ".tsv"))
# map_features = read_tsv("map_features_chr21.tsv")


#### get frequency of each trinucleotide per sequence (for 3rd process -- offset)

# add 1 bp downstream and upstream regions of width 1 or 2, so trinucs can be fetched
map_features = map_features %>%
  mutate(start = ifelse(width <= 2,
                        start-1,
                        start),
         end = ifelse(width <= 2,
                      end+1,
                      end),
         width = end-start+1)

# get sequence for each range
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                   names = makeGRangesFromDataFrame(map_features, 
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

# bind trinuc32 freqs to map_features
map_features_trinuc32_freq = map_features %>% 
  bind_cols(trinuc32_freq) %>% 
  group_by_at(vars(c(chromatin_features$name, AID_regions$name))) %>% 
  summarise_at(vars(matches("^[A,C,T,G][C,T][A,C,T,G]$")),
               ~sum(.)) %>% 
  mutate(chrom = chromosome) %>% 
  relocate(chrom)
gc()

write_tsv(map_features_trinuc32_freq, paste0("map_features_binarized_", chromosome, ".tsv"))
