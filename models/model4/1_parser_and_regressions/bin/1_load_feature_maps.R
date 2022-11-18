library(tidyverse)
# BiocManager::install('XVector', 'GenomicRanges', 'rtracklayer')
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")


### load chromatin feature maps and dna repair marks hg19 and keep as granges

# from command parameters
args = commandArgs(trailingOnly=TRUE)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.tsv",
                      no = args[1]) %>%
  read_tsv(comment = "#") 

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.tsv",
                            no = args[2]) %>%
  read_tsv(comment = "#")


dnarep_marks_list_files = list()
for (feature in dnarep_marks$name){
  
  path_file = filter(dnarep_marks, name == feature)$path
  
  dnarep_marks_list_files[[feature]] = tryCatch(import.bw(path_file),
                                                error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file),
                                                                                                      keep.extra.columns = T),
                                                                             error = function(e) import.bed(path_file)))
  # add feature name as the metadata (score) colname
  colnames(elementMetadata(dnarep_marks_list_files[[feature]])) = feature
  gc()
}

chromatin_features_list_files = list()
for (feature in chromatin_features$name){
  
  path_file = filter(chromatin_features, name == feature)$path
  
  chromatin_features_list_files[[feature]] = tryCatch(makeGRangesFromDataFrame(read_tsv(path_file),
                                                                               keep.extra.columns = T),
                                                      error = function(e) tryCatch(import.bed(path_file),
                                                                                   error = function(e) import.bw(path_file)))
  gc()
}


# merge Reptime and dna repair coordinates

dfleft = data.frame(chromatin_features_list_files[[1]]) %>%
  rename("chrom" = "seqnames")

list_feature_names = list()
list_feature_names[[1]] = colnames(elementMetadata(chromatin_features_list_files[[1]]))
n_chromatin_features = length(unlist(list_feature_names))

for (feature_i in seq(1, length(dnarep_marks_list_files))){
  
  dfright = data.frame(dnarep_marks_list_files[[feature_i]]) %>% rename("chrom" = "seqnames")
  
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

chr_names = paste0("chr", c(seq(1,22), "X", "Y"))
dfleft$chrom = factor(dfleft$chrom, ordered = T, levels = chr_names)

map_features = dfleft %>%
  rename("seqnames" = "chrom") %>%
  arrange(seqnames, start)

rm(dfleft) ; gc()

write_tsv(map_features, "map_features.tsv")
