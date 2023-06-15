library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")


args = commandArgs(trailingOnly=TRUE)

dnarep_mark = ifelse(interactive(),
                          yes = "H3K36me3",
                          no = args[1])

dnarep_mark_path = ifelse(interactive(),
                      	  yes = "/g/strcombio/fsupek_data/users/malvarez/git/trinucleotide_matching_tracks/matched_beds_symlinks/H3K36me3_pooled8_merged_2bins__3ntMatched_euclidean-0.1737.bed.gz",
                      	  no = args[2]) %>% 
  # remove comment
  gsub("( |\t).*", "", .)


# load file
dnarep_mark_file = tryCatch(import.bw(dnarep_mark_path),
                            error = function(e) tryCatch(import.bedGraph(dnarep_mark_path),
                                                         error = function(e) tryCatch(import.bed(dnarep_mark_path),
                                                                                      error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(dnarep_mark_path), keep.extra.columns = T),
                                                                                                                   error = function(e) import.wig(dnarep_mark_path)))))
# calculate median score for the dna repair mark of this iteration, for later binarization
scores = elementMetadata(dnarep_mark_file)[[1]]
if(is.numeric(scores)){
  median_score = median(scores)
  }else{
  median_score = paste(sort(unique(scores)), collapse = ",")
  }
  
write_tsv(data.frame(dnarep_mark, median_score),
          paste0("median_score_", dnarep_mark, ".tsv"))
