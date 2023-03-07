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
#counts = read_tsv("test_counts.tsv")
trinuc_sampling_per_bin_mymod_marina = function(counts, 
                                                stoppingCriterion = 0.001, # maximum tolerated distance (in relative frequency) in any column in any row; default 0.1% should be okay, if a bit stringent, for the 32 contexts (trinucleotide, strand-symmetrical)
                                                maxIter = 50000, # to prevent endless loops (in reality this can be a very high #, this works quite fast)
                                                n_finish_tokens = 1000){
  original_counts = counts
  tolerance = 10000000000000000000
  mintolerance = tolerance
  counts_mintolerance = counts
  finish_token = 0
  iter=0
  
  start_time = Sys.time()

  ### fix that counts offender row and correctable column are the ones intended, now it uses indices which is wrong since they are the indices of counts_not_too_few_mut_bins
  ### fix while loop that only ends when we have offender row with correctable column
  
  while ( TRUE ) {
    
    iter = iter + 1

    ## mark which bins have too few mutations, to ignore them for offender search (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
    too_few_mut_bins = filter(data.frame(rowSums(counts)),
                              `rowSums.counts.` <= (sum(rowSums(counts)) / length(rownames(counts)) / 100)) %>% # WARNING: div. by 100 so that not too many bins are excluded, revisit this
      rownames
    
    # check whether every row is declared as unusable (would result in an emtpy 'offender')
    if(length(too_few_mut_bins) >= length(rownames(counts))){
      counts = counts_mintolerance
      cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - 'too_few_mut_bins' comprises all possible bins, so no offender bin can be defined - Tolerance: %f -- Exiting and returning mintolerance results (%f)\nAnalysis terminated after %s\n", iter, maxIter, as.numeric(tolerance), mintolerance, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    # all rows are used for the "baseline mean frequencies"
    meanFreqs = colMeans(rowNorm(counts), na.rm = T)
    
    # ...but the 'too_few_mut_bins' are excluded from being candidates for offender,
    freqs_not_too_few_mut_bins = rownames_to_column(counts, "bin") %>% 
      ## so filter out bins (rows) with too few mutations
      filter(! bin %in% too_few_mut_bins) %>% 
      column_to_rownames("bin") %>% 
      rowNorm()
    
    got_offender_row_corr_col = F
    
    # LOOP HERE UNTIL WE HAVE AN OFFENDER ROW WITH CORRECTABLE COLUMN
    while (got_offender_row_corr_col = F){
      
      # find the 'offender' row, which is the row most different from the meanFreqs vector
      offender_name = which.max( apply(freqs_not_too_few_mut_bins, 1, function(x){ euclidean(x, meanFreqs) })) %>% # note that the meanFreqs DOES use ALL rows
        names()
      offender_row = rownames_to_column(counts, "bin") %>% 
        filter(bin == offender_name) %>% 
        column_to_rownames("bin")
      offender_row_freqs = rowNorm(offender_row)

      diffs = meanFreqs - offender_row_freqs
      
      if(sum(offender_row) <= 0){
        counts = counts_mintolerance
        cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - no more nts at the offender bin '%s' - Tolerance: %f -- Exiting and returning mintolerance results (%f)\nAnalysis terminated after %s\n", iter, maxIter, offender_name, as.numeric(tolerance), mintolerance, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        break
      }
  
      # in that row, find the column which is most responsible for the difference
      # however importantly we care ONLY about the negative differences in this vector!
      # i.e. those are the cases where the offending row has HIGHER freqs
      # (meaning we can correct that by removing sites... we can't add sites!!)
      
      is.correctableCol.zero = T
      
      while(is.correctableCol.zero == T){
        
        correctableCol = which.min(diffs)
  
        if(counts_not_too_few_mut_bins[offender, correctableCol] <= 0) {
          cat( sprintf("WARNING: At iter %i/%i, counts exhausted at bin '%s's 'correctableCol' (trinuc %s), trying next `which.min(diffs)` trinuc - Tolerance: %f\n%s have passed\n", iter, maxIter, offender, names(correctableCol), as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
          
          diffs[[correctableCol]] = NA
          
        } else {
          is.correctableCol.zero = F
          
          # check if we have removed all trinuc from diffs (so we should move on to next row to be used as offender)
          if(length(diffs) == 0){
  
            counts_not_too_few_mut_bins = c(counts_not_too_few_mut_bins, offender)
            
            cat(sprintf("WARNING: At iter %i/%i, removed all trinuc from 'diffs'; Bin '%s' is not used anymore as 'offender' bin - Tolerance: %f\n%s have passed\n", iter, maxIter, offender, as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
          } else {
            # end loop
            got_offender_row_corr_col = T
          }
        }
      }
  
      if(length(rownames(counts)) - length(too_few_mut_bins) <= 0){
        cat( sprintf("WARNING: At iter %i/%i, there are no bins left for being used as 'offender' bin! Resetting the 'counts_not_too_few_mut_bins' list so it can loop again from the bins with most mutations - Tolerance: %f\n%s have passed\n", iter, maxIter, as.numeric(tolerance), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        counts_not_too_few_mut_bins = c()
      }
    }

    worstCol = which.max(abs(diffs))  # this is sometimes the same as the correctable col
    
    # calculate tolerance -- note that tolerance is expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal
    tolerance = abs(diffs[worstCol]) 
    
    # store counts table if mintolerance is new minimum tolerance
    if(tolerance < mintolerance){
      mintolerance = tolerance
      counts_mintolerance = counts
    }
    
    # did we reduce the difference enough?
    if (tolerance <= stoppingCriterion) {
      cat(sprintf("Successfully completed optimization: tolerance reached. Poorest match at row %d col %d - Returning current results\n", offender, worstCol) );
      break
    }
    
    # note this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
    subtractThis = round( diffs[correctableCol] * sum(counts[offender, ]) )
    
    if ( subtractThis == 0 ) {
      # try ×'n_finish_tokens' times more
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
    
    if (iter == maxIter) {
      counts = counts_mintolerance
      cat( sprintf("Stopping optimization - maximum number of iterations reached - Returning mintolerance results (%f)\n", mintolerance) );
      break
    }  
    
    ## output log every 100th iter
    if(iter %% 100 == 0){
      cat( sprintf("Running iter %i/%i... Subtracted %i muts from trinuc %s at bin '%s' - Euclidean: %.3f, Tolerance: %f, ShortestWin: %7d, AvgWin: %7d\n%s have passed\n", iter, maxIter, abs(as.numeric(subtractThis)), names(correctableCol), offender, euclidean(freqs[offender,],meanFreqs), as.numeric(tolerance), min( rowSums(counts) ), round(mean( rowSums(counts) )), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
    }

  } # keep iterating...
  
  return(rownames_to_column(counts, "bin"))
}
