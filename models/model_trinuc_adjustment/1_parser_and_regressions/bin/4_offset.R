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
                                yes = lapply(list(c(
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr1.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr2.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr3.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr4.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr5.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr6.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr7.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr8.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr9.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr10.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr11.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr12.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr13.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr14.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr15.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr16.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr17.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr18.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr19.tsv"))[1],
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr20.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr21.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr22.tsv"))[1])),
                                  #Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chrX.tsv"))[1])),
                                  read_tsv),
                                no = lapply(list(args), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) 

gc()


### calculate offset (log(n trinuc of each of the 32 types (e.g. ACT) that exist in each RT-dnarepmarks combination, and could therefore be any A(C>D)T SNV))

offset_temp = map_features_binarized %>%
  # sum up trinuc32 frequencies within repliseq-dnamarks profile
  select(-chrom) %>% 
  # RepliSeq == 0 has only NNNNNN.. sequences, so no trinucs are found.. so remove Repliseq==0 bin
  filter(RepliSeq != 0) %>% 
  group_by_at(vars(!matches("^[A,C,G,T][C,T][A,C,G,T]$"))) %>% 
  summarise_at(vars(matches("^[A,C,G,T][C,T][A,C,G,T]$")),
               ~ sum(.)) %>% 
  ungroup %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "freq_trinuc32") %>% 
  # triplicate each row, adding '( >A)', '( >G)' and '( >T)' around the central C or T
  group_by_at(vars(!matches("^freq_trinuc32$"))) %>% 
  slice(rep(row_number(), 3))

AGT_column = rep(c('>A)', '>G)', '>T)'), 
                 times = length(rownames(offset_temp)) / 3) %>% 
  data.frame %>% 
  `colnames<-`("AGT")

offset_temp = offset_temp %>% 
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
  relocate(tri, .before = freq_trinuc32)

## add an offset==0 for all RT×dnamarks combinations that do not exist, e.g. if in `RT_1 × ogg30_high × ogg60_low × UV1_high × UV2_high × MSH6_low × ...` genome regions there are no SNVs, offset==0
# NOTE: usually this is not needed, as all possible combinations exist, but in that case the merge(offset, empty_offset, all = T) does nothing
ncols = length(names(offset_temp)) - 3 # 3 == tri + freq_trinuc32 + repliseq columns
row_high = c("1", rep("high", ncols), "A(C>A)A")
row_low = c("2", rep("low", ncols), "A(C>G)A")
col_repliseq = data.frame(seq(2, 6, 1)) %>% `colnames<-`("RepliSeq") %>% mutate(RepliSeq = as.character(RepliSeq))
col_tri = data.frame(unique(offset_temp$tri)) %>% `colnames<-`("tri") %>% filter(!tri %in% c("A(C>A)A", "A(C>G)A"))

empty_offset = data.frame(matrix(ncol = length(names(offset_temp)))) %>% 
  `colnames<-`(names(offset_temp)) %>% 
  select(-c(freq_trinuc32)) %>% 
  rbind(row_high) %>% 
  rbind(row_low) %>% 
  drop_na %>% 
  merge(col_repliseq, all = T) %>% 
  merge(col_tri, all = T) %>% 
  expand(!!! syms(names(offset_temp)[!str_detect(names(offset_temp), "freq_trinuc32")])) %>%
  drop_na

if(!is.null(empty_offset$AID_regions)){
  empty_offset = empty_offset %>% 
    mutate(AID_regions = gsub("low", "AID_target", AID_regions),
           AID_regions = gsub("high", "bgGenome", AID_regions))
}

offset = merge(offset_temp, empty_offset, all = T) %>% 
  replace_na(list(freq_trinuc32 = 0)) %>% 
  mutate(trinuc32 = gsub("\\(", "", tri),
         trinuc32 = gsub("\\>.*\\)", "", trinuc32)) %>% 
  relocate(trinuc32, .before = tri)

write_tsv(offset, "offset.tsv")
