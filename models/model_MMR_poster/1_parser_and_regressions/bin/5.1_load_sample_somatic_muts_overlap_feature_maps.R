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


### load sample; merge with dna repair, chromatin landscape, parse

# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample = ifelse(interactive(),
                yes = "6faf3507-bacc-5782-b834-f48c89fc10c7", #"MSM0.103", #"MSM0.124",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/parse_vcfs/res/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/downsample_mutburden/res/muts_pass_",
                                no = args[2]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

dnarep_marks = ifelse(interactive(),
                      yes = "../input_lists/dnarep_marks.csv",
                      no = args[3]) %>%
  read_csv(comment = "#")

chromatin_features = ifelse(interactive(),
                            yes = "../input_lists/chromatin_features.csv",
                            no = args[4]) %>%
  read_csv(comment = "#")


## load map_features (SINGLE chromosome) from 2nd process
dfleft = ifelse(interactive(),
                yes = Sys.glob("map_features_chr21.tsv")[1],
                no = args[5]) %>% 
  fread %>% as_tibble %>% 
  rename("chrom" = "seqnames")
gc()

chromosome = unique(dfleft$chrom)


## NEW keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                 yes = "/g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed",
                                 no = args[6]) %>%
  import.bed() %>% data.frame %>%
  mutate(seqnames = gsub("^", "chr", seqnames)) %>%
  rename("chrom" = "seqnames") %>% 
  filter(chrom == chromosome)


# load collected median_scores from 1st process
median_scores = ifelse(interactive(),
                       yes = lapply(list(c(Sys.glob("median_score_DHS.tsv")[1],
                                           Sys.glob("median_score_exons.tsv")[1],
                                           Sys.glob("median_score_H3K36me3.tsv")[1],
                                           Sys.glob("median_score_RepliSeq.tsv")[1],
                                           Sys.glob("median_score_RnaSeq.tsv")[1])), 
                                    read_tsv),
                       no = lapply(list(args[-(1:6)]), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .)


## load sample (somatic mutations)

# find path+sample name that exists
existing_file = c()
for(file_exists in paste0(path_somatic_variation, sample, ".csv")){
  if(file.exists(file_exists)){
    existing_file = c(existing_file, file_exists)
  }
}

# raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
if(length(existing_file) != 1){
  stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
}

# load it
dfright = read_csv(existing_file) %>%
  select(chr, start, end, tri) %>% 
  rename("chrom" = "chr") %>%
  filter(chrom == chromosome) %>% 
  mutate(mut_id = paste0("mut_", row_number())) %>% 
  # NEW keep SNVs in good mappability regions
  bed_intersect(good_mappability_regions, suffix = c("_dfright", "_crg75")) %>% 
  select(c("chrom", contains("_dfright"))) %>% 
  rename_all(~str_replace_all(., "_dfleft|_dfright", ""))

gc()

### map chromatin features
merged = dfright %>%
  # intersect with dfleft
  bed_intersect(dfleft, suffix = c("_dfright", "_dfleft")) %>%
  select(-c(contains("chrom"), contains("start_"), contains("end_"), contains("width_"), contains("strand_"))) %>%
  # remove the "dleft" and "dright" parts of the column names
  rename_all(~str_replace_all(., "_dfleft|_dfright", "")) %>%
  # combine chromatin features, if there is any (although there should typically be only RepliSeq)
  {if (length(chromatin_features$name)>0) unite("mb_domain", contains(chromatin_features$name), sep = "_") else .} %>%
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
  mutate(sample_id = sample)
gc()

if(nrow(merged) == 0){
  sprintf("WARNING - Empty 'merged' table: probably 'good_mappability_regions' has removed all SNVs for this sample! Table will be empty...")
}

write_tsv(merged, paste0("ready_for_regression_", sample, "_", chromosome, ".tsv"))
