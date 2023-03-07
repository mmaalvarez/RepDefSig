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
                          yes = "TP53_dauno_K562",
                          no = args[1])

dnarep_mark_path = ifelse(interactive(),
                      	  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/TP53/gaps_score0/GSM3781507_WT_KO_Dauno_coverage_gaps-to-score-0.bed",
                      	  no = args[2]) %>% 
  # remove comment
  gsub("( |\t).*", "", .)


# load file
dnarep_mark_file = tryCatch(import.bw(dnarep_mark_path),
                            error = function(e) tryCatch(import.bedGraph(dnarep_mark_path),
                                                         error = function(e) tryCatch(import.bed(dnarep_mark_path),
                                                                                      error = function(e) read_tsv(dnarep_mark_path))))

# calculate median score for the dna repair mark of this iteration, for later binarization
scores = elementMetadata(dnarep_mark_file)[[1]]
if(is.numeric(scores)){
  median_score = median(scores)
  }else{
  median_score = paste(sort(unique(scores)), collapse = ",")
  }
  

write_tsv(data.frame(dnarep_mark, median_score),
          paste0("median_score_", dnarep_mark, ".tsv"))
