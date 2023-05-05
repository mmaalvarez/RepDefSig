library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
library(rlang)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")


# - For each DNA repair pathway (e.g. MMR):
# 	- 1st) calculate the mutation burden baseline (=="baseline sample"): mean # mutations (treating each of the 96 trinuc types separately, to keep relative %s thereof) across all control samples, WITHIN each of the 2 genome BINS (i.e. high/low MMR repair abundance)
# 	- 2nd) each simulated sample will be the "baseline sample" + an increased # mutations (proportional to mutfoldinc), but ONLY in the genome bin with 'high' activity of a mark for that repair (e.g. high MSH6)
# 		- again, each trinuc type treated independently to keep the relative trinuc %s
#     - UPDATE: mutfoldinc multiplied to the baseline sample's genome-wide mutation burden per trinuc (e.g. if baseline sample's total # mutations at ATA in target bin (e.g. AID-targets) is just 2, and in bg bin (bgGenome) is 250, then multiply to the target bin a FRACTION of the total # ATA muts of the bg bin: e.g. a 0.01 --> 2 * 0.01 * 250
# - All the simulated samples will be then run together with real ones in the NMF/VAE
# - The simulated samples should have high exposures of the signatures with high weights for the expected mark (e.g. for "MMRdef simulated samples", high "MSH6-weight signature exposure")
# 	- the more sensitive the model, low mut. burden-increase (e.g. ×2?) simulated samples should already show high exposures


# from command parameters
args = commandArgs(trailingOnly=TRUE)

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_",
                                no = args[1]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

# controls for creating the 'baseline sample', for then adding to it extra mutations to create simulated positive controls
control_samples = ifelse(interactive(),
                         yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv",
                         no = args[2]) %>%
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  map_df(~read_tsv(.x)) %>% 
  filter(str_detect(info2, "[C,c]ontrol")) %>% 
  pull(sample_id)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[3]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[4]) %>%
  read_csv(comment = "#")

# load offset from 4th process
offset = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1],
                no = args[5]) %>% 
  read_tsv
# rename the chromatin environment column (typically 'RepliSeq') to match the "mb_domain" name given to the general mutation table
colnames(offset)[1] = "mb_domain"
offset = offset %>% 
  mutate(mb_domain = factor(mb_domain))

## load map_features (SINGLE chromosome) from 2nd process
dfleft = ifelse(interactive(),
                yes = Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/map_features_chr1.tsv"),
                no = args[6]) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()

chromosome = unique(dfleft$chrom)


## NEW keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                 yes = "/g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed",
                                 no = args[7]) %>% 
  import.bed() %>% data.frame %>%
  mutate(seqnames = gsub("^", "chr", seqnames)) %>% 
  rename("chrom" = "seqnames") %>% 
  filter(chrom == chromosome)


## load "mutfoldinc"
mutfoldinc = ifelse(interactive(),
                    yes = "0.01",
                    no = args[8]) %>% 
  as.numeric

# which dna repair mark has muts increased in this iteration
dnarep_mark_simulate = ifelse(interactive(),
                               yes = "OGG1_GOx30_chipseq",
                               no = args[9])

## load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c(Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_OGG1_GOx30_chipseq.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_OGG1_GOx60_chipseq.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_UV_XRseq_NHF1_PP64_1h_Rep1.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_UV_XRseq_NHF1_CPD_1h.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_XRCC4.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_SETD2_control.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_MSH6_control.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_TP53_dauno_K562.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_TP53_dauno_MOLM13.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_AID_regions.tsv")[1],
                                           Sys.glob("../work/[[:alnum:]][[:alnum:]]/*/median_score_A3A_TpCpH_hairpins.tsv")[1])), 
                                    read_tsv),
                       no = lapply(list(args[-(1:9)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)



## load control samples' somatic mutations

for(control_sample in control_samples){
  
  # find path+sample name that exists
  existing_file = c()
  for(file_exists in paste0(path_somatic_variation, control_sample, ".csv")){
    if(file.exists(file_exists)){
      existing_file = c(existing_file, file_exists)
    }
  }
  
  # raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
  if(length(existing_file) != 1){
    stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
  }
  
  ## load sample
  dfright = read_csv(existing_file) %>%
    select(chr, start, end, tri) %>% 
    rename("chrom" = "chr")  %>% 
    filter(chrom == chromosome) %>% 
    mutate(mut_id = paste0("mut_", row_number())) %>%
    # NEW keep SNVs in good mappability regions
    bed_intersect(good_mappability_regions, suffix = c("_dfright", "_crg75")) %>% 
    select(c("chrom", contains("_dfright"))) %>% 
    rename_all(~str_replace_all(., "_dfleft|_dfright", ""))
  
  ### if there are mutations at this chromosome, map chromatin features
  if(length(rownames(dfright)) >= 1){

    ## map chromatin features
    merged = dfright %>%
      # intersect with dfleft
      bed_intersect(dfleft, suffix = c("_dfright", "_dfleft")) %>%
      select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
      # remove the "dleft" and "dright" parts of the column names
      rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
      # combine chromatin features (although there should typically be only RepliSeq)
      unite("mb_domain", contains(chromatin_features$name), sep = "_") %>% 
      # binarize weighted average DNA repair value by being lower or larger than the across-genome median
      rowwise %>% 
      mutate_at(vars(contains(match = dnarep_marks$name)),
                function(x){var_name = rlang::as_label(substitute(x))
                ifelse(!is.na(suppressWarnings(as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score))),
                       # it's numeric score, binarize as low/high
                       yes = ifelse(x <= as.numeric(filter(median_scores, dnarep_mark == var_name)$median_score),
                                    yes = "low",
                                    no = "high"),
                       # it's factor, leave as is
                       no = x)}) %>% 
      ### mut count
      select(-mut_id) %>% 
      table %>%
      as.data.frame %>%
      rename("mutcount" = "Freq") %>% 
      # add offset (table) to have 0 muts in the non-existing combinations...
      merge(offset, all = T) %>%
      replace_na(list(mutcount = 0)) %>% 
      relocate(mutcount) %>% 
      # ... but not really using the actual offset column at this point
      select(-freq_trinuc32) %>% 
      mutate("control_sample" = control_sample)
    
    if(nrow(merged) == 0){
      stop("ERROR - Empty 'merged' table: probably 'good_mappability_regions' has removed all SNVs for this CONTROL sample! Exiting...\n")
    }
    
    ### append to control_samples_table if this exists, otherwise create it
    if(exists("control_samples_table")){
      control_samples_table = bind_rows(control_samples_table, merged)
    }else{
      control_samples_table = merged
    }
    rm(merged)
    
  } else {
    
    cat(sprintf("WARNING: Sample %s does not have mutations at %s, so it will not be used for creating the 'baseline sample' for this chromosome...\n", control_sample, chromosome))
  }
  
  gc()
}


### create baseline sample

# first exclude control outliers (too many mutations, as maybe there are cryptic mutagens)
mut_burden_controls = control_samples_table %>% 
  group_by(control_sample) %>% 
  summarise(mutburden = sum(mutcount)) %>% 
  pull(mutburden) %>% 
  quantile

if((mut_burden_controls[["75%"]] + 1.5 * (mut_burden_controls[["75%"]] - mut_burden_controls[["25%"]]))  >  max(mut_burden_controls)){
  cat(sprintf("There are no 'extreme' outliers\n"))
}

control_samples_no_outliers = control_samples_table %>% 
  group_by(control_sample) %>% 
  summarise(mutburden = sum(mutcount)) %>% 
  # use as cutoff for removing treated samples the mean of mut_burden_controls
  filter(mutburden <= mean(mut_burden_controls)) %>% 
  pull(control_sample)

baseline_sample = control_samples_table %>% 
  filter(control_sample %in% control_samples_no_outliers) %>% 
  select(-c(control_sample, .overlap)) %>% 
  group_by(across(!contains("mutcount"))) %>% 
  summarise(mean_mutcount = mean(mutcount)) %>% 
  ungroup %>% 
  mutate(chr = chromosome,
         `mutfoldinc` = mutfoldinc,
         simulated_mark = dnarep_mark_simulate)

write_tsv(baseline_sample, paste0("baseline_sample_", chromosome, "_mutfoldinc", mutfoldinc, "_", dnarep_mark_simulate, ".tsv"))
