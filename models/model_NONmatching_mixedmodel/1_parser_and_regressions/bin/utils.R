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


##########################################################################################################################
## function for parsing all track files: load, rename chromosomes, extract iteration's chromosome, and add +-50bp to ranges if there are gaps
# imported by 2_load_feature_maps.R

parse_feature_files = function(features_table){
  
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
  
    ## NEW if the ranges are not continuous, i.e. there are gaps...
    range_gaps = gaps(feature_file) %>% 
      data.frame %>% 
      filter(strand == "*" & seqnames == chromosome & start >= 2)
    
    if(nrow(range_gaps) >= 1){
  
      cat(sprintf("%s contains gaps; adding 50bp on both sides of its ranges...\n", feature))
      
      ## ...shorten gaps between ranges (+-50bp) to ensure more SNVs overlap with this repmark
      
      # (BUT IF the adjacent range is CLOSER than 100bp, just join them)
      if(min(data.frame(range_gaps)$width) <= 101){
        
        cat(sprintf("Some of %s's gaps are shorter than 101bp, this will take longer...\n", feature))
        
        # functions
        width_upstream_gap = function(range_start_pos){
          filter(range_gaps, end == range_start_pos-1) %>% 
            pull(width)}
        width_downstream_gap = function(range_end_pos){
          filter(range_gaps, start == range_end_pos+1) %>% 
            pull(width)}
        
        feature_file = data.frame(feature_file) %>% 
          arrange(start) %>% 
          rowwise() %>% 
          mutate(start = ifelse(length(width_upstream_gap(start)) != 0,
                                # there is a gap ending before "start" of this range
                                yes = ifelse(width_upstream_gap(start) >= 101,
                                             # the gap's width is >=100bp, there is room for extending the range 50bp without overlapping the previous range
                                             yes = start - 50,
                                             # the gap's width is <100bp, so only extend half of its width
                                             no = start - floor(width_upstream_gap(start) / 2)),
                                # there is no gap ending before "start" (e.g. no gap between ranges, or start of chromosomes?)
                                no = start),
                 # now do the equivalent thing for 'end'
                 end = ifelse(length(width_downstream_gap(end)) != 0,
                              # there is a gap starting after the "end" of this range
                              yes = ifelse(width_downstream_gap(end) >= 101,
                                           # the gap's width is >=100bp, there is room for extending the range 50bp without overlapping the next range
                                           yes = end + 50,
                                           # the gap's width is <100bp, so only extend half of its width
                                           no = end + floor(width_downstream_gap(end) / 2)),
                              # there is no gap starting after the "end" (e.g. no gap between ranges, or end of chromosomes?)
                              no = end)) %>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
        
      } else { # every gap >100bp --> just add +-50bp to every range
        
        cat(sprintf("None of %s's gaps are shorter than 101bp, this will be fast!\n", feature))
        
        feature_file = data.frame(feature_file) %>% 
          arrange(start) %>% 
          mutate(start = start - 50,
                 end = end + 50) %>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
      }
    } else {
      cat(sprintf("%s does not contain gaps, add to list and continue...\n", feature))
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
# imported by 3_binarize_scores.R

trinuc_sampling_per_bin = function(counts, 
                                   stoppingCriterionTolerance = 0.001,  # this is the stopping cutoff for the "tolerance" (i.e. the highest (positive) trinuc frequency difference between current offender bin vs. "target" bin (average of all other bins) -- e.g. GCT being 21.7% of target's trinucs and 2.9% of offender's could be a tolerance of 0.188); it might be too low i.e. too many mutations lost to reach it (when subsampling from introns it is feasible, but in my case maybe not)
                                   stoppingCriterionVar_score = 0.1,    # NEW this is for the var_score, see below
                                   min_rm_muts = 100){
  
  # initialize functions/parameters/lists/constants
  stoppingCriterionTolerance = as.numeric(stoppingCriterionTolerance)
  stoppingCriterionVar_score = as.numeric(stoppingCriterionVar_score)
  speed_adjustment = 1 # NEW to speed up adjustment, see below -- starts as 1 but then gets increased
  min_rm_muts = as.numeric(min_rm_muts) # NEW also to speedup the adjustment, i.e. this is the min. num of mutations removed
  if(stoppingCriterionVar_score < 1){
    excess_var_cutoff = as.numeric(gsub("1e-", "", as.character(format(stoppingCriterionVar_score, scientific = T)))) + 1 # NEW in case var_score is excess_var_cutoff times larger than the registered minVarScore (lowest var_score), quit and return the output with that minVarScore -- excess_var_cutoff is the number of decimals in stoppingCriterionVar_score
  } else {
    stop("ERROR: 'stoppingCriterionVar_score' has to be <1! Exiting...\n")
  }
  maxIter = 40000*length(counts)  # to prevent endless loops (in reality this can be a very high #, this works quite fast)
  trinuc_fraction_zero_muts = 0.5 # see below
  bin_with_freqs_like_rest = c() # to store bin names whose trinuc freqs DIFFERENCES with respect to the freqs of all the remaining bins combined are already too low, to not use them anymore as offender bin
  iter = 1
  proportion_of_total_muts_per_bin = rowSums(counts) / sum(rowSums(counts)) # NEW this is to don't remove too many mutations from the bins with more mutations, because if a bin has originally e.g. 30% of total mutations we don't want that after subsampling it has e.g. 10% of total mutations, since that wouldn't be a good representation of the effect of the bin's properties (dna repair abundances, RT) on mutation burden

  # function to normalize each row to sum to 1
  rowNorm = function(m) {
    t( apply(m, 1, function(x) { x/sum(x) } ) )
  }
  
  start_time = Sys.time() # track time
  
  while ( TRUE ) {

    ## mark which bins have too few mutations, to ignore them for offender search (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
    too_few_mut_bins = filter(data.frame(rowSums(counts)),
                              `rowSums.counts.` <= (sum(rowSums(counts)) / length(rownames(counts)))) %>% 
      rownames
    
    if(length(unique(c(bin_with_freqs_like_rest, too_few_mut_bins))) >= length(rownames(counts))){
      # this declares every row as unusable (would result in an emtpy 'offender'), so reinitiate the least critical of the 2 restrictions
      bin_with_freqs_like_rest = c()
      cat(sprintf("WARNING: At iter %i/%i, 'bin_with_freqs_like_rest' and 'too_few_mut_bins' comprise all possible bins, so no offender bin can be defined; Restarting 'bin_with_freqs_like_rest' - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
    }
    
    ## setting the 'offender' row name to the bin with MOST mutations summing all 32 trinuc
    offender = rowSums(counts) %>% 
      data.frame %>% 
      arrange(desc(`.`)) %>% 
      ## make sure bins already in same freq that remaining bins (bin_with_freqs_like_rest) are not considered
      rownames_to_column("bin") %>% 
      filter(! bin %in% bin_with_freqs_like_rest) %>% 
      ## ALSO filter out bins (rows) with too few mutations (i.e. less than the total current sum of mutations (across rows AND columns) divided by the number of bins)
      filter(! bin %in% too_few_mut_bins) %>% 
      column_to_rownames("bin") %>% 
      slice_head(n = 1) %>% 
      rownames()
    if(sum(counts[offender,]) <= 0){
      counts = counts_minVarScore
      cat(sprintf("ERROR: Cannot continue optimization at iter %i/%i - no more nts at the offender bin '%s' - Tolerance: %f, Var score: %f -- Exiting and returning the table with minVarScore (%f)\nAnalysis terminated after %s\n", iter, maxIter, offender, as.numeric(tolerance), var_score, minVarScore, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }
    
    ## setting the 'target' table to the sum of ALL the remaining bins, since instead of just 2 rows (exons and introns) we have many (1 per genomic bin) 
    target = filter(counts, row.names(counts) != offender) %>% 
      summarise_all(sum)
    
    ## compute the normalized frequencies & compare them between the 'target bin' vs the offender bin
    freqs_target = rowNorm(target)
    
    freqs_all_rows = data.frame(rowNorm(counts))
    freqs_offender = freqs_all_rows[offender,]
    
    diffs = freqs_target - freqs_offender
    
    
    ## in that bin (row), find the column which is most responsible for the difference; however importantly we care ONLY about the negative differences in this vector! i.e. those are the cases where the offending row has HIGHER freqs (meaning we can correct that by removing sites... we can't add sites!!)
    diffs_correctableCol = diffs
    
    is.correctableCol.zero = T
    
    while(is.correctableCol.zero == T){
      
      correctableCol = which.min(diffs_correctableCol)

      if(counts[offender, correctableCol] <= 0) {
        cat( sprintf("WARNING: At iter %i/%i, counts exhausted at bin '%s's 'correctableCol' (trinuc %s), trying next `which.min(diffs)` trinuc - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, offender, names(correctableCol), as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
        
        diffs_correctableCol[[correctableCol]] = NULL
        
      } else {
        is.correctableCol.zero = F
        
        # check if we have removed all trinuc from diffs_correctableCol (so we should move on to next row to be used as offender)
        if(length(diffs_correctableCol) == 0){

          bin_with_freqs_like_rest = c(bin_with_freqs_like_rest, offender)
          
          cat(sprintf("WARNING: At iter %i/%i, removed all trinuc from 'diffs_correctableCol'; Added bin '%s' to bin_with_freqs_like_rest, so it's not used anymore as 'offender' bin - # bins not in this list: %d - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, offender, length(rownames(counts)) - length(bin_with_freqs_like_rest), as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        }
      }
    }

    
    ### computing tolerance -- expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal
    worstCol = which.max(abs(diffs))  # this is sometimes the same as the correctable col
    worstCol_tri = colnames(counts)[worstCol]
    tolerance = abs(diffs[worstCol])
    if(iter == 1){initTolerance = tolerance}
    
    
    ### NEW also computing 'var_score': an additional measure of mean (across trinucs) mut freq variances across rows (bins), it has to be minimised as well
    var_score = freqs_all_rows %>% 
      # multiply freqs (sum up to 1 per bin) by 32trinuc (i.e. so it sums up to 32 per bin)
      mutate_all(~.*32) %>% 
      summarise_all(~var(., na.rm = T)) %>%  # na.rm in case a bin (row) had all 0s before rowNorm, remove it since the NaNs that result from it dont allow to summarise
      as.numeric %>% mean
    
    if(iter == 1){
      # in 1st iteration, store initial var_score and initialize min var score
      initVarScore = var_score
      minVarScore = initVarScore
      
    } else { # remaining iterations

      # check if we have reached the lowest var_score so far
      if(var_score < minVarScore){
        # update min score
        minVarScore = var_score
        # store the current counts table since it has reached the lowest var_score yet
        counts_minVarScore = counts
      }      
      
      # did we reduce the tolerance AND var_score enough already?
      if (tolerance <= stoppingCriterionTolerance & var_score <= stoppingCriterionVar_score) {
        cat(sprintf("Successfully completed optimization: tolerance and var_score cutoffs reached at iter %i/%i, when %i muts were subtracted from trinuc %s at bin '%s' - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, abs(as.numeric(subtractThis)), names(correctableCol), offender, as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        break
      }
      
      # if the var_score is excess_var_cutoff times larger than the minVarScore, quit and return the last counts_minVarScore instead of the current counts table
      if (var_score >= 1+minVarScore*excess_var_cutoff) {
        
        # first try a last time by doing speed_adjustment==2
        if(speed_adjustment<2){
          speed_adjustment = 2
          # also allow double excess_var_cutoff
          excess_var_cutoff = excess_var_cutoff*2
          cat(sprintf("WARNING: At iter %i/%i, var_score > minVarScore*%f (%i muts were subtracted from trinuc %s at bin '%s' - Tolerance: %f, Var score: %f) -- Trying with speed_adjustment==2\n%s have passed\n", iter, maxIter, excess_var_cutoff, abs(as.numeric(subtractThis)), names(correctableCol), offender, as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        } else {
          # if speed_adjustment==2 already tried, exit
          counts = counts_minVarScore
          cat(sprintf("WARNING: At iter %i/%i, var_score > minVarScore*%f (%i muts were subtracted from trinuc %s at bin '%s' - Tolerance: %f, Var score: %f) -- Exiting and returning the table with minVarScore (%f)\n%s have passed\n", iter, maxIter, excess_var_cutoff, abs(as.numeric(subtractThis)), names(correctableCol), offender, as.numeric(tolerance), var_score, minVarScore, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
          break
        }
      }
    }
    
    ### this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
    adjustment = (diffs[correctableCol] * sum(counts[offender, ]) + 1) * speed_adjustment
    subtractThis = round(adjustment - adjustment*proportion_of_total_muts_per_bin[[offender]])
    
    ## NEW instead of `cat( sprintf("Successfully completed optimization: count reduction <=0.5. Poorest match at bin %s trinuc %s (%d)\n", offender, worstCol_tri, diffs[[worstCol_tri]]) ); break;`
    if ( subtractThis == 0 ) {
      # store the name of this bin to not use it anymore as offender
      bin_with_freqs_like_rest = c(bin_with_freqs_like_rest, offender)
      
      # to remove 1 mut, not 0...
      subtractThis = 1
      
      if(length(bin_with_freqs_like_rest) %% 10 == 0){
        cat( sprintf("WARNING: At iter %i/%i, added bin '%s' to 'bin_with_freqs_like_rest', so it's not used anymore as 'offender' bin - # bins not in this list: %d - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, offender, length(rownames(counts)) - length(bin_with_freqs_like_rest), as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      }
      if(length(rownames(counts)) - length(too_few_mut_bins) - length(bin_with_freqs_like_rest) <= 0){
        cat( sprintf("WARNING: At iter %i/%i, there are no bins left for being used as 'offender' bin! Resetting the 'bin_with_freqs_like_rest' list so it can loop again from the bins with most mutations - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
        bin_with_freqs_like_rest = c()
      }
    }

    ##### decrease counts in the responsible column to get closer to the rowNorm; 'min_rm_muts' is to speedup
    subtractThis = subtractThis - min_rm_muts
    counts[offender, correctableCol] = counts[offender, correctableCol] + subtractThis

    # just in case this resulted in negative mut counts, make them 0
    if(counts[offender, correctableCol] <= -1){
      counts[offender, correctableCol] = 0
    }
    
    # update and check number of iterations
    iter=iter+1
    if (iter==maxIter) {
      counts = counts_minVarScore
      cat( sprintf("WARNING: Stopping optimization - maximum number of iterations (%i) reached - Tolerance: %f, Var score: %f -- Exiting and returning the table with minVarScore (%f)\nAnalysis terminated after %s\n", maxIter, as.numeric(tolerance), var_score, minVarScore, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
      break
    }  
    ## output log every 100th iter
    if(iter %% 100 == 0){
      cat( sprintf("Running iter %i/%i... Subtracted %i muts from trinuc %s at bin '%s' - Tolerance: %f, Var score: %f\n%s have passed\n", iter, maxIter, abs(as.numeric(subtractThis)), names(correctableCol), offender, as.numeric(tolerance), var_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
    }
  }
  counts = rownames_to_column(counts, "bin")
  return(counts)
}


## test count table
# counts = Sys.glob("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/AID_SHM/OLD_no_flanks_offtargets/1_extract_SNVs_lymphoid_tumors__AID_on_and_off_target_regions_vs_whole_genome/res/*/output/SBS/*_AID2_SNVs.SBS96.all") %>%
#   lapply(read_tsv) %>%
#   Reduce(function(x, y) left_join(x, y, by="MutationType"), .) %>%
#   mutate(MutationType = gsub("\\[", "", MutationType),
#          MutationType = gsub(">.\\]", "", MutationType)) %>%
#   group_by(MutationType) %>%
#   summarise_all(.funs = sum) %>%
#   ungroup() %>%
#   pivot_longer(cols = c(contains("_On_off"), contains("_Unclustered")), names_to = "bin", values_to = "SNVs") %>%
#   mutate(SNVs = ifelse(str_detect(SNVs, "_Unclustered"),
#                        SNVs/1000,
#                        SNVs)) %>% 
#   pivot_wider(names_from = MutationType, values_from = SNVs) %>%
#   column_to_rownames("bin")

## dummy test count table
# counts = read_tsv("dummy.tsv") %>%
#   column_to_rownames("bin")
