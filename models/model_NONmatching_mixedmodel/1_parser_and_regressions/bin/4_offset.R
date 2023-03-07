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
                                yes = lapply(c(
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr1.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr2.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr3.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr4.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr5.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr6.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr7.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr8.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr9.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr10.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr11.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr12.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr13.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr14.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr15.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr16.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr17.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr18.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr19.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr20.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr21.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr22.tsv")),
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chrX.tsv"))),
                                  read_tsv),
                                no = lapply(args, read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) 

gc()


### calculate offset (log(n trinuc of each of the 32 types (e.g. ACT) that exist in each RT-dnarepmarks combination, and could therefore be any A(C>D)T SNV))

offset = map_features_binarized %>%
  # sum up trinuc32 frequencies within repliseq-dnamarks profile; add a dummy count; and log
  select(-chrom) %>% 
  group_by_at(vars(!matches("^[A,C,G,T][C,T][A,C,G,T]$"))) %>% 
  summarise_at(vars(matches("^[A,C,G,T][C,T][A,C,G,T]$")),
               ~ log(sum(.) + 1)) %>% 
  ungroup %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "log_freq_trinuc32") %>% 
  # RepliSeq == 0 has only NNNNNN.. sequences, so no trinucs are found.. so remove Repliseq==0 bin
  filter(RepliSeq != 0) %>% 
  # triplicate each row, adding '( >A)', '( >G)' and '( >T)' around the central C or T
  group_by_at(vars(!matches("^log_freq_trinuc32$"))) %>% 
  slice(rep(row_number(), 3))

AGT_column = rep(c('>A)', '>G)', '>T)'), 
                 times = length(rownames(offset)) / 3) %>% 
  data.frame %>% 
  `colnames<-`("AGT")

offset = offset %>% 
  bind_cols(AGT_column) %>% 
  mutate(tri = paste0(substr(trinuc32, start = 1, stop = 1),
                      "(",
                      substr(trinuc32, start = 2, stop = 2),
                      AGT,
                      substr(trinuc32, start = 3, stop = 3)),
         # correct T>T to T>C
         tri = gsub("T>T", "T>C", tri)) %>% 
  ungroup %>% 
  select(-c(trinuc32, AGT)) %>%
  relocate(tri, .before = log_freq_trinuc32)

## add an offset==0 for all RT×dnamarks×32tri combinations that do not exist, e.g. if in `RT_1 × ogg30_high × ogg60_low × UV1_high × UV2_high × MSH6_low × ...` genome regions there is no ACG, offset==0 because it's impossible that there is an A(C>D)G SNV in those regions
# NOTE: usually this is not needed, as all possible combinations exist, but in that case the merge(offset, empty_offset, all = T) does nothing
ncols = length(names(offset)) - 3 # 3 == log_freq_trinuc32 + tri + repliseq columns
row_high = c("1", rep("high", ncols), "A(C>A)A")
row_low = c("2", rep("low", ncols), "A(C>G)A")
col_repliseq = data.frame(seq(2, 6, 1)) %>% `colnames<-`("RepliSeq") %>% mutate(RepliSeq = as.character(RepliSeq))
col_tri = data.frame(unique(offset$tri)) %>% `colnames<-`("tri") %>% filter(!tri %in% c("A(C>A)A", "A(C>G)A"))

empty_offset = data.frame(matrix(ncol = length(names(offset)))) %>% 
  `colnames<-`(names(offset)) %>% 
  select(-log_freq_trinuc32) %>% 
  rbind(row_high) %>% 
  rbind(row_low) %>% 
  drop_na %>% 
  merge(col_repliseq, all = T) %>% 
  merge(col_tri, all = T) %>% 
  relocate(tri, .after = last_col()) %>% 
  expand(!!! syms(names(offset)[!str_detect(names(offset), "log_freq_trinuc32")])) %>%
  drop_na %>% 
  mutate(AID_regions = gsub("high", unique(offset$AID_regions)[1], AID_regions),
         AID_regions = gsub("low", unique(offset$AID_regions)[2], AID_regions))

offset = merge(offset, empty_offset, all = T) %>% 
  replace_na(list(log_freq_trinuc32 = 0))

# for 9 binary dna marks × 6 RT bins × 32 NC|TN trinuc32 == XXX --> × 3 AGT(C) == XXX rows
write_tsv(offset, "offset.tsv")
