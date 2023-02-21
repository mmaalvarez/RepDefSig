library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging (bed_intersect)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")

source("utils.R")


### load chromatin feature maps and dna repair marks hg19 and keep as granges

# from command parameters
args = commandArgs(trailingOnly=TRUE)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[1]) %>%
  read_csv(comment = "#") 

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[2]) %>%
  read_csv(comment = "#")

chromosome = ifelse(interactive(),
                    yes = "chr21",
                    no = paste0("chr", args[3]))


### apply function from utils.R (for parsing all track files: load, rename chromosomes, extract iteration's chromosome, and add +-50bp to ranges if there are gaps)

dnarep_marks_list_files = parse_feature_files(dnarep_marks) # %>% filter(name %in% c('XRCC4', 'AID_regions')))

chromatin_features_list_files = parse_feature_files(chromatin_features)



### merge Reptime and dna repair coordinates

if(length(chromatin_features_list_files) == 1){
  dfleft = data.frame(chromatin_features_list_files[[1]]) %>% rename("chrom" = "seqnames")
}else{
  exit("ERROR: There is not just 1 chromatin feature (e.g. RepliSeq); there are either none or >1! Exiting...\n")
}

list_feature_names = list()
list_feature_names[[1]] = colnames(elementMetadata(chromatin_features_list_files[[1]]))
n_chromatin_features = length(unlist(list_feature_names))

for (feature_i in seq(1, length(dnarep_marks_list_files))){

  dfright = data.frame(dnarep_marks_list_files[[feature_i]]) %>% 
    rename("chrom" = "seqnames")
  
  list_feature_names[[n_chromatin_features + feature_i]] = colnames(elementMetadata(dnarep_marks_list_files[[feature_i]]))

  dfleft = bed_intersect(dfleft, dfright, suffix = c("_dfleft",
                                                     "_dfright"))
  gc()
  rm(dfright)
  gc()
  dfleft = dfleft %>% 
    mutate(start = ifelse(start_dfleft >= start_dfright,
                          start_dfleft,
                          start_dfright),
           end = ifelse(end_dfleft <= end_dfright,
                        end_dfleft,
                        end_dfright),
           strand = "*",
           width = end - start + 1) %>%
    select(chrom, start, end, width, strand, contains(unlist(list_feature_names)))
  gc()

  colnames(dfleft) = gsub("_dfleft", "", colnames(dfleft))
  colnames(dfleft) = gsub("_dfright", "", colnames(dfleft))
  gc()
}
rm(dnarep_marks_list_files)
gc()

dfleft = dfleft %>%
  rename("seqnames" = "chrom") %>%
  arrange(start)
gc()

# this will go to next (per chromosome) and last (all chromosomes collected) processes
write_tsv(dfleft, paste0("map_features_", chromosome, ".tsv"))
