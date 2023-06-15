library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging (bed_intersect)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("desc", "dplyr")
conflict_prefer("reverseComplement", "spgs")

#setwd("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin")


##########################################################################################################################
## function for parsing all track files: load, rename chromosomes, extract iteration's chromosome, and POSSIBLY add +- 'nt_extend' bp to ranges if there are gaps
# imported by 2_load_feature_maps.R

parse_feature_files = function(features_table, 
                               chromosome,
                               nt_extend = as.numeric(0)){
  ## auxiliar functions
  width_upstream_gap = function(range_start_pos){
    filter(range_gaps, end == range_start_pos-1) %>% 
      pull(width)}

  width_downstream_gap = function(range_end_pos){
    filter(range_gaps, start == range_end_pos+1) %>% 
      pull(width)}

  if(nt_extend >= 1){
    cat(sprintf("MODE: Removing %fbp from both sides of range gaps\n", nt_extend))
  } else {
    cat(sprintf("MODE: NOT removing bp from both sides of range gaps\n"))
  }
  
  output_list = list()
  
  for (feature in features_table$name){
  
    path_file = filter(features_table, name == feature)$path
  
    feature_file = tryCatch(import.bw(path_file),
                            error = function(e) tryCatch(import.bedGraph(path_file),
                                                         error = function(e) tryCatch(import.bed(path_file),
                                                                                      error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file), keep.extra.columns = T),
                                                                                                                  error = function(e) import.wig(path_file)))))
    # keep only metadata without 0s
    for(metadata_col in names(elementMetadata(feature_file))){
      if(!"FALSE" %in% unique(unique(elementMetadata(feature_file)[[metadata_col]]) == 0)){
        mcols(feature_file)[[metadata_col]] <- NULL
      }
    }
    # make sure that there is only 1 (and not more) metadata columns
    if(length(names(elementMetadata(feature_file))) != 1){
      stop(paste0("Feature file ", feature, " has 0 or >1 metadata columns! Exiting...\n"))
    }

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
    
    feature_file = feature_file[seqnames(feature_file) == chromosome]
    
    # this extra step is to remove + and - strands from some features such as UV ones, as they are not used
    feature_file = feature_file %>% 
      data.frame %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
  
    ## ONLY do all this if we wanted to extend range extremes, to reduce gaps between them
    if(nt_extend >= 1){
      
      ## if the ranges are not continuous, i.e. there are gaps...
      range_gaps = gaps(feature_file) %>% 
        data.frame %>% 
        filter(strand == "*" & seqnames == chromosome & start >= 2)
      
      if(nrow(range_gaps) >= 1){
    
        cat(sprintf("%s contains gaps; adding %fbp on both sides of its ranges...\n", feature, nt_extend))
        
        ## ...shorten gaps between ranges (+- nt_extend bp) to ensure more SNVs overlap with this repmark
        
        # (BUT IF the adjacent range is CLOSER than nt_extend*2 + 1 bp, just join them)
        if(min(data.frame(range_gaps)$width) <= nt_extend*2 + 1){
          
          cat(sprintf("Some of %s's gaps are shorter than %fbp, this will take longer...\n", feature, nt_extend*2+1))
          
          feature_file = data.frame(feature_file) %>% 
            arrange(start) %>% 
            rowwise() %>% 
            mutate(start = ifelse(length(width_upstream_gap(start)) != 0,
                                  # there is a gap ending before "start" of this range
                                  yes = ifelse(width_upstream_gap(start) >= nt_extend*2 +1,
                                               # the gap's width is >= nt_extend*2 +1 bp, there is room for extending the range nt_extend bp without overlapping the previous range
                                               yes = start - nt_extend,
                                               # the gap's width is < nt_extend*2 +1 bp, so only extend half of its width
                                               no = start - floor(width_upstream_gap(start) / 2)),
                                  # there is no gap ending before "start" (e.g. no gap between ranges, or start of chromosomes?)
                                  no = start),
                   # now do the equivalent thing for 'end'
                   end = ifelse(length(width_downstream_gap(end)) != 0,
                                # there is a gap starting after the "end" of this range
                                yes = ifelse(width_downstream_gap(end) >= nt_extend*2 +1,
                                             # the gap's width is >= nt_extend*2 +1 bp, there is room for extending the range nt_extend bp without overlapping the next range
                                             yes = end + nt_extend,
                                             # the gap's width is < nt_extend*2 +1 bp, so only extend half of its width
                                             no = end + floor(width_downstream_gap(end) / 2)),
                                # there is no gap starting after the "end" (e.g. no gap between ranges, or end of chromosomes?)
                                no = end)) %>%
            makeGRangesFromDataFrame(keep.extra.columns = T)
          
          } else { # every gap > nt_extend*2+1 bp --> just add +- nt_extend bp to every range
            
            cat(sprintf("None of %s's gaps are shorter than %fbp, this will be fast!\n", feature, nt_extend*2 +1))
          
            feature_file = data.frame(feature_file) %>% 
              arrange(start) %>% 
              mutate(start = start - nt_extend,
                     end = end + nt_extend) %>%
              makeGRangesFromDataFrame(keep.extra.columns = T)
            }
        } else {
          cat(sprintf("%s does not contain gaps, add to list and continue...\n", feature))
        }
    }
    
    output_list[[feature]] = feature_file
    rm(feature_file)
    gc()
  }
  return(output_list)
}
