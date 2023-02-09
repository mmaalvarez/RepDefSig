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
                          yes = "OGG1_GOx30_chipseq",
                          no = args[1])

dnarep_mark_path = ifelse(interactive(),
                      	  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/DNA_repair/BER/OGG1/chipseq/bed_hg19/hg19_GSM2357435_CP-Sample_Flag-OGG1-GOX-30-1.tdf.bed #             BER/OGG1/GOx30        chip-seq  .bed",
                      	  no = args[2]) %>% 
  # remove comment
  gsub("( |\t).*", "", .)


low_mappability_regions = read_tsv(args[3], col_names = F)


# load file
dnarep_mark_file = tryCatch(import.bw(dnarep_mark_path),
                            error = function(e) tryCatch(import.bedGraph(dnarep_mark_path),
                                                         error = function(e) tryCatch(read_tsv(dnarep_mark_path),
                                                                                      error = function(e) import.bed(dnarep_mark_path))))

# calculate median score for the dna repair mark of this iteration, for later binarization
median_score = median(elementMetadata(dnarep_mark_file)[[1]])

write_tsv(data.frame(dnarep_mark, median_score),
          paste0("median_score_", dnarep_mark, ".tsv"))
