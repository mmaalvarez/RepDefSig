library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
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
conflict_prefer("expand", "tidyr")


# load map_features_binarized (all chromosomes) from 2nd process
args = commandArgs(trailingOnly=TRUE)

map_features_binarized = ifelse(interactive(),
                                yes = lapply("map_features_binarized_chr21.tsv", read_tsv),
                                no = lapply(args, read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) 

gc()


### calculate offset (log(n trinuc of each of the 32 types (e.g. ACT) that exist in each RT-dnarepmarks combination, and could therefore be any A(C>D)T SNV))

offset = map_features_binarized %>%
  # sum up ACROSS ALL MATCHED trinuc32 frequencies within repliseq-dnamarks profile; add a dummy count; and log
  rename("chrom" = "seqnames") %>% 
  select(-chrom) %>% 
  mutate(log_freq_trinuc32 = log(rowSums(select(., matches("^[A,C,G,T][C,T][A,C,G,T]$"))) + 1)) %>% 
  # RepliSeq == 0 has only NNNNNN.. sequences, so no trinucs are found.. so remove Repliseq==0 bin
  filter(RepliSeq != 0) %>% 
  select(-matches("^[A,C,G,T][C,T][A,C,G,T]$"))

## add an offset==0 for all RT×dnamarks combinations that do not exist, e.g. if in `RT_1 × ogg30_high × ogg60_low × UV1_high × UV2_high × MSH6_low × ...` genome regions there are no SNVs, offset==0
# NOTE: usually this is not needed, as all possible combinations exist, but in that case the merge(offset, empty_offset, all = T) does nothing
ncols = length(names(offset)) - 2 # 2 == log_freq_trinuc32 + repliseq columns
row_high = c("1", rep("high", ncols))
row_low = c("2", rep("low", ncols))
col_repliseq = data.frame(seq(2, 6, 1)) %>% `colnames<-`("RepliSeq") %>% mutate(RepliSeq = as.character(RepliSeq))

empty_offset = data.frame(matrix(ncol = length(names(offset)))) %>% 
  `colnames<-`(names(offset)) %>% 
  select(-log_freq_trinuc32) %>% 
  rbind(row_high) %>% 
  rbind(row_low) %>% 
  drop_na %>% 
  merge(col_repliseq, all = T) %>% 
  expand(!!! syms(names(offset)[!str_detect(names(offset), "log_freq_trinuc32")])) %>%
  drop_na

if(!is.null(empty_offset$AID_regions)){
  empty_offset = empty_offset %>% 
    mutate(AID_regions = gsub("low", "AID_target", AID_regions),
           AID_regions = gsub("high", "bgGenome", AID_regions))
}

offset = merge(offset, empty_offset, all = T) %>% 
  replace_na(list(log_freq_trinuc32 = 0)) %>% 
  rename("mb_domain" = "RepliSeq")

write_tsv(offset, "offset.tsv")
