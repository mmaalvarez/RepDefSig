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
conflict_prefer("desc", "dplyr")


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
    cat(sprintf("MODE: Removing %ibp from both sides of range gaps\n", nt_extend))
  } else {
    cat(sprintf("MODE: NOT removing bp from both sides of range gaps\n"))
  }
  
  output_list = list()
  
  for (feature in features_table$name){
  
    path_file = filter(features_table, name == feature)$path
  
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
    
        cat(sprintf("%s contains gaps; adding %ibp on both sides of its ranges...\n", feature, nt_extend))
        
        ## ...shorten gaps between ranges (+- nt_extend bp) to ensure more SNVs overlap with this repmark
        
        # (BUT IF the adjacent range is CLOSER than nt_extend*2 + 1 bp, just join them)
        if(min(data.frame(range_gaps)$width) <= nt_extend*2 + 1){
          
          cat(sprintf("Some of %s's gaps are shorter than %ibp, this will take longer...\n", feature, nt_extend*2+1))
          
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
            
            cat(sprintf("None of %s's gaps are shorter than %ibp, this will be fast!\n", feature, nt_extend*2 +1))
          
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




###########################################################################################################################

### function for matching mutated trinucleotide (32) proportions between genomic bins
# takes mutation count dataframe with N rows (1 per genomic bin, e.g. RT6-OGG1low-UVhigh-...) × 32 trinuc type columns
# removes mut counts, so that mutated trinucleotide proportions are "matched" across bins
# it also removes accordingly the total trinuc counts that exist in that bin, so that the offset will also account for the matching
# modified from marina/fran's
# can be imported by 5.2 and 6.2

trinuc_matching = function(total_trinuc_table, mut_trinuc_table, 
                           stoppingCriterion = 0.01, # desired Euclidean score (max. distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies)
                           maxIter = 20000*length(total_trinuc_table), # to prevent endless loops (in reality this can be a very high #, this works quite fast)
                           n_finish_tokens = 1000,
                           max_fraction_removed_muts = 0.25){
  ## auxiliar functions
  rowNorm = function(m){t(apply(m, 1, function(x){x/sum(x)}))} # normalize each row to sum to 1
  euclidean = function(a, b){sqrt(sum((a-b)^2))} 
  
  ## initialize constants/variables
  euclidean_score = 10000000000000000000 # init as an absurdly large value
  mineuclidean_score = euclidean_score
  counts = list("total_trinuc" = total_trinuc_table, 
                "mut_trinuc" = mut_trinuc_table)
  counts_mineuclidean = counts
  finish_token = 0
  removed_muts = 0
  total_orig_muts = sum(rowSums(counts$mut_trinuc))
  max_removed_muts = total_orig_muts * max_fraction_removed_muts
  iter = 0
  start_time = Sys.time()

  while ( TRUE ) {
    
    iter = iter + 1
    
    ## check if we have reached max. num. iterations
    if (iter == maxIter) {
      counts = counts_mineuclidean
      cat( sprintf("Stopping optimization - maximum number of iterations reached - Returning min Euclidean score results (%f)\n", mineuclidean_score) );
      break
    }  
    
    ## check if we have already removed too many mutations
    if(removed_muts > max_removed_muts){
      counts = counts_mineuclidean
      cat(sprintf("Stopping optimization at iter %i/%i - %.02f%% of the %i original mutations have already been removed -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, removed_muts/total_orig_muts*100, total_orig_muts, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    ## mark which bins have too few mutations, to ignore them for offender search (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
    too_few_mut_bins = filter(data.frame(rowSums(counts$mut_trinuc)),
                              `rowSums.counts.mut_trinuc.` <= (sum(rowSums(counts$mut_trinuc)) / length(rownames(counts$mut_trinuc)) / 2)) %>% # WARNING: div. by 10 so that not too many bins are excluded, revisit this
      rownames
    
    ## check whether every row has been declared as unusable (would result in an emtpy 'offender')
    if(length(rownames(counts$mut_trinuc)) - length(too_few_mut_bins) <= 0){
      counts = counts_mineuclidean
      cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - 'too_few_mut_bins' comprises all possible bins, so no offender bin can be defined -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    ###################################################
    ## get OFFENDER ROW WITH CORRECTABLE COLUMN
    
    ## all rows are used for the "baseline mean frequencies"
    meanFreqs = colMeans(rowNorm(counts$mut_trinuc), na.rm = T)
    # ...but the 'too_few_mut_bins' are excluded from being candidates for offender...
    freqs_not_too_few_mut_bins = rownames_to_column(counts$mut_trinuc, "bin") %>% 
      ## ...so filter out bins (rows) with too few mutations
      filter(! bin %in% too_few_mut_bins) %>% 
      column_to_rownames("bin") %>% 
      rowNorm()
    
    ## LOOP HERE UNTIL WE HAVE AN OFFENDER ROW WITH CORRECTABLE COLUMN
    
    got_offender_row_corr_col = F
    
    while (got_offender_row_corr_col == F){
      
      # find the 'offender' row, which is the row most different from the meanFreqs vector
      offender_name = which.max(apply(freqs_not_too_few_mut_bins, 1, function(x){ euclidean(x, meanFreqs) })) %>% # note that the meanFreqs DOES use ALL rows
        names()
      offender_row = rownames_to_column(counts$mut_trinuc, "bin") %>% 
        filter(bin == offender_name) %>% 
        column_to_rownames("bin")
      offender_row_freqs = rowNorm(offender_row)

      diffs = data.frame(meanFreqs - offender_row_freqs)
      
      if(sum(offender_row) <= 0){
        counts = counts_mineuclidean
        cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - no more nts at the offender bin '%s' -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, offender_name, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        break
      }
  
      ## in that row, find the column which is most responsible for the difference; however importantly we care ONLY about the negative differences in this vector! i.e. those are the cases where the offending row has HIGHER freqs (meaning we can correct that by removing sites... we can't add sites!!)
      
      is.correctableCol.zero = T
      
      while(is.correctableCol.zero == T){
        
        correctableColIndex = which.min(diffs)
        correctableCol = names(correctableColIndex)
        correctableColCounts = select(offender_row, all_of(correctableCol))
  
        if(correctableColCounts <= 0) {
          cat( sprintf("WARNING: At iter %i/%i, counts exhausted at bin '%s's 'correctableCol' (trinuc %s), trying next `which.min(diffs)` trinuc - Euclidean score: %f\n%s have passed\n", iter, maxIter, offender_name, correctableCol, as.numeric(euclidean_score), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
          
          diffs = mutate_at(diffs, 
                            vars(all_of(correctableCol)), 
                            ~gsub(".*", NA, .))
          
        } else {
          is.correctableCol.zero = F
          
          # check if we have removed all trinuc from diffs (so we should move on to next row to be used as offender)
          if(length(diffs) == 0){
  
            # don't use this bin anymore as offender
            too_few_mut_bins = c(too_few_mut_bins, offender_name)
            
            cat(sprintf("WARNING: At iter %i/%i, removed all trinuc from 'diffs'; Bin '%s' is not used anymore as 'offender' bin - Euclidean score: %f\n%s have passed\n", iter, maxIter, offender_name, as.numeric(euclidean_score), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
          
          } else {
            # end loop
            got_offender_row_corr_col = T
          }
        }
      }
  
      if(length(rownames(counts$mut_trinuc)) - length(too_few_mut_bins) <= 0){
        counts = counts_mineuclidean
        cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - 'too_few_mut_bins' comprises all possible bins, so no offender bin can be defined -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        break
      }
    }
    
    ## calculate euclidean score (how diff. are the offender row freqs. from the mean freqs. across the table)
    euclidean_score = euclidean(offender_row_freqs, meanFreqs)
    
    # store counts table if euclidean_score is new minimum
    if(euclidean_score < mineuclidean_score){
      mineuclidean_score = euclidean_score
      counts_mineuclidean = counts
    }
    
    # did we reduce the difference enough?
    if (euclidean_score <= stoppingCriterion) {
      cat(sprintf("Successfully completed optimization: Euclidean score (%f) lower than %f - Returning current results\n", euclidean_score, stoppingCriterion) )
      break
    }
  

    ###############################
    ### subtract mutations from a trinuc of the offender bin
    
    ### note this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
    subtractThis = as.numeric(round(diffs[correctableColIndex] * sum(offender_row)))
    
    if ( subtractThis == 0 ) {
      # try ×'n_finish_tokens' times more
      subtractThis = 10
      finish_token = finish_token + 1
      cat( sprintf("Count reduction <= 0.5 - Euclidean score: %f - %i/%i tokens used, so continue...\n", euclidean_score, finish_token, n_finish_tokens) )
      
      if(finish_token >= n_finish_tokens){
        counts = counts_mineuclidean
        cat( sprintf("Stopping optimization - count reduction <= 0.5 and all %i tokens have been used - Returning min Euclidean score results (%f)\n", n_finish_tokens, mineuclidean_score) )
        break
      }
    }
    
    ## now simply decrease counts in the responsible column to get closer to the mean
    counts$mut_trinuc = rownames_to_column(counts$mut_trinuc, "bin") %>% 
      mutate_at(vars(all_of(correctableCol)),
                ~ifelse(bin == offender_name,
                        . + subtractThis,
                        .)) %>% 
      column_to_rownames("bin")
    
    # also at the corresponding total trinuc counts
    counts$total_trinuc = rownames_to_column(counts$total_trinuc, "bin") %>% 
      mutate_at(vars(all_of(correctableCol)),
                ~ifelse(bin == offender_name,
                        . + subtractThis,
                        .)) %>% 
      column_to_rownames("bin")
    
    removed_muts = removed_muts + abs(subtractThis)
    

    ## output log every 1000th iter
    if(iter %% 1000 == 0){
      cat( sprintf("Running iter %i/%i... Subtracted %i muts from trinuc %s at bin '%s' - %.02f%% muts have been removed - Euclidean score: %f\n%s have passed\n", iter, maxIter, abs(as.numeric(subtractThis)), correctableCol, offender_name, removed_muts/total_orig_muts*100, euclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
    }
    
    ## keep iterating...
  }
  
  ## return final counts_mineuclidean tables (both the total trinuc and the mutated trinuc tables)
  return(lapply(counts, 
                rownames_to_column, "bin"))
}



###############################################################################
## ALTERNATIVE to MATCHING: adjustment based on dividing # mutated trinucleotides by the total (# mutated + # non-mutated + 1 dummy count) per bin, sum fractions, and multiply by the original total # mutations in that bin
# faster (no iterations)
# should keep signal (ratios of mutations between bins) less altered
# can be imported by 5.2 and 6.2

trinuc_adjustment = function(total_trinuc_table, mut_trinuc_table){
  
  adjusted_mutcounts = round(rowSums(mut_trinuc_table/(total_trinuc_table+1) * rowSums(mut_trinuc_table))) %>% 
    data.frame %>% 
    `colnames<-`("mutcount") %>% 
    rownames_to_column("bin")
  
  return(adjusted_mutcounts)
}
