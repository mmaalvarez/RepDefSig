library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(rlang)
library(conflicted)
conflict_prefer("slice", "dplyr")


dnarep_marks_orig = read_csv("../dnarep_marks_orig_paths.csv", comment = "#") %>% 
  filter(name != "AID_regions")

dnarep_marks_list = list()

for (feature in dnarep_marks_orig$name){
  
  path_file = filter(dnarep_marks_orig, name == feature)$path
  
  feature_file = tryCatch(import.bw(path_file),
                          error = function(e) tryCatch(import.bedGraph(path_file),
                                                       error = function(e) tryCatch(import.bed(path_file),
                                                                                    error = function(e) makeGRangesFromDataFrame(read_tsv(path_file),
                                                                                                                                 keep.extra.columns = T))))
  
  # add feature name as the metadata (score) colname
  colnames(elementMetadata(feature_file)) = feature
  
  # extract chromosome (args[3]), otherwise it gets too long and the bed_intersect crashes
  chr_list = levels(seqnames(feature_file))
  if(!str_detect(paste(chr_list, collapse="|"),"chr[1-9]")){
    cat(sprintf("Feature '%s' chromosome names lack a '^chr' prefix; adding it...\n", feature))
    feature_file = data.frame(feature_file) %>% 
      mutate(seqnames = gsub("^", "chr", seqnames)) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
  }
  
  # this extra step is to remove + and - strands from some features such as UV ones, as they are not used
  feature_file = feature_file %>% 
    data.frame %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)

  dnarep_marks_list[[feature]] = feature_file
  rm(feature_file)
  gc()
}

set.seed(1)

for(feature in dnarep_marks_orig$name){
  
  dnarep_marks_list[[feature]] = data.frame(dnarep_marks_list[[feature]]) %>% 
    mutate(!!feature := sample(get(feature))) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  gc()
}

for(feature in dnarep_marks_orig$name){
  
  gr = dnarep_marks_list[[feature]]
  df = data.frame(seqnames = seqnames(gr),
                  start = start(gr)-1,
                  end = end(gr),
                  score = mcols(gr))
  
  write_tsv(df, paste0(feature, "_shuffled_scores.bed"), col_names = F)
}
