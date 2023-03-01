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

# helper functions
width_upstream_gap = function(range_start_pos){
  filter(range_gaps, end == range_start_pos-1) %>% 
    pull(width)}

width_downstream_gap = function(range_end_pos){
  filter(range_gaps, start == range_end_pos+1) %>% 
    pull(width)}

# main function
parse_feature_files = function(features_table, 
                               chromosome,
                               nt_extend = as.numeric(0)){
  
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
### trinucleotide (32) sampling per genomic bin
# takes mutation count dataframe with N rows (1 per genomic bin, e.g. RT6-OGG1low-UVhigh-...) × 32 trinuc type columns
# returns the downsampled frequencies for each trinuc
# 2 versions: trinuc_sampling_per_bin_mymod_marina() and trinuc_sampling_per_bin_mymod_ahmed()
# imported by 3_binarize_scores.R


## helper functions for trinuc_sampling_per_bin_mymod_marina

# jitter should be between 0 (no jitter) and 1 (100% jitter)
generateJitteryRow = function(freqs, numNt, jitter) {
  trunc( freqs * runif(length(freqs), min=1-jitter, max=1+jitter) * numNt)
}

# normalize each row to sum to 1. There may be a faster way.
rowNorm =function(m) {
  t( apply(m, 1, function(x) { x/sum(x) } ) )
}

euclidean = function(a, b) {
  sqrt(sum((a - b)^2))
} 

## main function (GOOD, new)
trinuc_sampling_per_bin_mymod_marina = function(counts, 
                                                stoppingCriterion = 0.001, # maximum tolerated distance (in relative frequency) in any column in any row; default 0.1% should be okay, if a bit stringent, for the 32 contexts (trinucleotide, strand-symmetrical)
                                                maxIter = 50000, # to prevent endless loops (in reality this can be a very high #, this works quite fast)
                                                n_finish_tokens = 1000){
  original_counts = counts
  mintolerance = 10000000000000000000
  counts_mintolerance = counts
  finish_token = 0
  iter=0
  
  start_time = Sys.time()

  while ( TRUE ) {
    
    iter=iter+1;

    ## mark which bins have too few mutations, to ignore them for offender search (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
    too_few_mut_bins = filter(data.frame(rowSums(counts)),
                              `rowSums.counts.` <= (sum(rowSums(counts)) / length(rownames(counts)) / 1000)) %>% 
      rownames
    
    # check whether every row is declared as unusable (would result in an emtpy 'offender')
    if(length(too_few_mut_bins) >= length(rownames(counts))){
      counts = counts_mintolerance
      cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - 'too_few_mut_bins' comprises all possible bins, so no offender bin can be defined - Tolerance: %f -- Exiting and returning mintolerance results (%f)\nAnalysis terminated after %s\n", iter, maxIter, as.numeric(tolerance), mintolerance, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    ## filter out bins (rows) with too few mutations (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
    counts_not_too_few_mut_bins = rownames_to_column(counts, "bin") %>% 
      filter(! bin %in% too_few_mut_bins) %>% 
      column_to_rownames("bin")

    freqs=rowNorm(counts_not_too_few_mut_bins);
    meanFreqs=colMeans(freqs, na.rm = T);
    
    # find the 'offending' row which is most different from the mean relFreq vector
    offender = which.max( apply(freqs, 1, function(x){ euclidean(x,meanFreqs) })  );

    if(sum(counts[offender,]) <= 0){
      counts = counts_mintolerance
      cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - no more nts at the offender bin '%s' - Tolerance: %f -- Exiting and returning mintolerance results (%f)\nAnalysis terminated after %s\n", iter, maxIter, offender, as.numeric(tolerance), mintolerance, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    diffs = meanFreqs - freqs[offender,];
    
    # in that row, find the column which is most responsible for the difference
    # however importantly we care ONLY about the negative differences in this vector!
    # i.e. those are the cases where the offending row has HIGHER freqs
    # (meaning we can correct that by removing sites... we can't add sites!!)
    
    is.correctableCol.zero = T
    
    while(is.correctableCol.zero == T){
      
      correctableCol = which.min(diffs)

      if(counts[offender, correctableCol] <= 0) {
        cat( sprintf("WARNING: At iter %i/%i, counts exhausted at bin '%s's 'correctableCol' (trinuc %s), trying next `which.min(diffs)` trinuc - Tolerance: %f\n%s have passed\n", iter, maxIter, offender, names(correctableCol), as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
        
        diffs[[correctableCol]] = NA
        
      } else {
        is.correctableCol.zero = F
        
        # check if we have removed all trinuc from diffs (so we should move on to next row to be used as offender)
        if(length(diffs) == 0){

          counts_not_too_few_mut_bins = c(counts_not_too_few_mut_bins, offender)
          
          cat(sprintf("WARNING: At iter %i/%i, removed all trinuc from 'diffs'; Bin '%s' is not used anymore as 'offender' bin - Tolerance: %f\n%s have passed\n", iter, maxIter, offender, as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        }
      }
    }

    if(length(rownames(counts)) - length(too_few_mut_bins) <= 0){
      cat( sprintf("WARNING: At iter %i/%i, there are no bins left for being used as 'offender' bin! Resetting the 'counts_not_too_few_mut_bins' list so it can loop again from the bins with most mutations - Tolerance: %f\n%s have passed\n", iter, maxIter, as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      counts_not_too_few_mut_bins = c()
    }

    worstCol = which.max(abs(diffs))  # this is sometimes the same as the correctable col
    
    # calculate tolerance -- note that tolerance is expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal
    tolerance = abs(diffs[worstCol]) 
    
    # store counts table if mintolerance
    if(tolerance < mintolerance){
      mintolerance = tolerance
      counts_mintolerance = counts
    }
    
    # did we reduce the difference enough?
    if ( tolerance <= stoppingCriterion ) {
      cat( sprintf("Successfully completed optimization: tolerance reached. Poorest match at row %d col %d - Returning current results\n", offender, worstCol) );
      break;
    }
    
    # note this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
    subtractThis = round( diffs[correctableCol] * sum(counts[offender, ]) )
    
    if ( subtractThis == 0 ) {
      # try n times more
      subtractThis = 10
      finish_token = finish_token+1
      cat( sprintf("Count reduction <=0.5. Poorest match at row %d col %d - Continue (finish_token nº%i)...\n", offender, worstCol, finish_token) )
      
      if(finish_token >= n_finish_tokens){
        counts = counts_mintolerance
        cat( sprintf("Successfully completed optimization: count reduction <=0.5. Poorest match at row %d col %d - Returning mintolerance results (%f)\n", offender, worstCol, mintolerance) )
        break
      }
    }
    
    # now simply decrease counts in the responsible column to get closer to the mean
    counts[offender, correctableCol] = counts[offender, correctableCol] + subtractThis
    
    if (iter==maxIter) {
      counts = counts_mintolerance
      cat( sprintf("Stopping optimization - maximum number of iterations reached - Returning mintolerance results (%f)\n", mintolerance) );
      break;
    }  
    
    ## output log every 100th iter
    if(iter %% 100 == 0){
      cat( sprintf("Running iter %i/%i... Subtracted %i muts from trinuc %s at bin '%s' - Euclidean: %.3f, Tolerance: %f, ShortestWin: %7d, AvgWin: %7d\n%s have passed\n", iter, maxIter, abs(as.numeric(subtractThis)), names(correctableCol), offender, euclidean(freqs[offender,],meanFreqs), as.numeric(tolerance), min( rowSums(counts) ), round(mean( rowSums(counts) )), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
    }

  } # keep iterating...
  
  return(rownames_to_column(counts, "bin"))
}
