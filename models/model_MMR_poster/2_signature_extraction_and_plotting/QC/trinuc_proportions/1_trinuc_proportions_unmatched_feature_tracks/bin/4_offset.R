library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg19") # all coordinates are in hg19
library(rlang)
library(spgs)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")


# load map_features_binarized (all chromosomes) from 2nd process
args = commandArgs(trailingOnly=TRUE)

map_features_binarized = ifelse(interactive(),
                                yes = lapply(list(c(
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr1.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr2.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr3.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr4.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr5.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr6.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr7.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr8.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr9.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr10.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr11.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr12.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr13.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr14.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr15.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr16.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr17.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr18.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr19.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr20.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr21.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chr22.tsv"))[1],
                                  Sys.glob(paste0("../work/[[:alnum:]][[:alnum:]]/*/map_features_binarized_chrX.tsv"))[1])),
                                  read_tsv),
                                no = lapply(list(args), read_tsv)) %>%
  Reduce(function(x, y) bind_rows(x, y), .) 
gc()


### calculate offset (log(n trinuc of each of the 32 types (e.g. ACT) that exist in each RT-dnarepmarks combination, and could therefore be any A(C>D)T SNV))
offset_temp = map_features_binarized %>%
  # sum up trinuc32 frequencies within repliseq-dnamarks profile
  select(-chrom) %>% 
  # RepliSeq == 0 has only NNNNNN.. sequences, so no trinucs are found.. so remove Repliseq==0 bin (IF RepliSeq EXISTS)
  filter(if (exists("RepliSeq", where = cur_data())) RepliSeq != 0 else TRUE) %>% 
  group_by_at(vars(!matches("^[A,C,G,T][C,T][A,C,G,T]$"))) %>% 
  summarise_at(vars(matches("^[A,C,G,T][C,T][A,C,G,T]$")),
               ~ sum(.)) %>% 
  ungroup %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "freq_trinuc32") %>% 
  # triplicate each row, adding '( >A)', '( >G)' and '( >T)' around the central C or T
  group_by_at(vars(!matches("^freq_trinuc32$"))) %>% 
  slice(rep(row_number(), 3))

AGT_column = rep(c('>A)', '>G)', '>T)'), 
                 times = length(rownames(offset_temp)) / 3) %>% 
  data.frame %>% 
  `colnames<-`("AGT")

offset_temp = offset_temp %>% 
  bind_cols(AGT_column) %>% 
  mutate(tri = paste0(substr(trinuc32, start = 1, stop = 1),
                      "(",
                      substr(trinuc32, start = 2, stop = 2),
                      AGT,
                      substr(trinuc32, start = 3, stop = 3)),
         # correct T>T to T>C
         tri = gsub("T>T", "T>C", tri)) %>% 
  ungroup %>% 
  select(-c(trinuc32, AGT)) %>%
  relocate(tri, .before = freq_trinuc32)

## add an offset==0 for all RT×dnamarks combinations that do not exist, e.g. if in `RT_1 × ogg30_high × ogg60_low × UV1_high × UV2_high × MSH6_low × ...` genome regions there are no SNVs, offset==0
# NOTE: usually this is not needed, as all possible combinations exist, but in that case the merge(offset, empty_offset, all = T) does nothing

if(("RepliSeq" %in% colnames(offset_temp))  &  (0 %in% unique(offset_temp$RepliSeq))){
  ncols = length(names(offset_temp)) - 3 # 3 == tri + freq_trinuc32 + Repliseq columns
  row_high = c("1", rep("high", ncols), "A(C>A)A")
  row_low = c("2", rep("low", ncols), "A(C>G)A")
  col_tri = data.frame(unique(offset_temp$tri)) %>% `colnames<-`("tri") %>% filter(!tri %in% c("A(C>A)A", "A(C>G)A"))
  col_repliseq = data.frame(seq(2, 6, 1)) %>% `colnames<-`("RepliSeq") %>% mutate(RepliSeq = as.character(RepliSeq))
  
} else {
  # no RepliSeq (or there is RepliSeq, but not with a "0" level)
  ncols = length(names(offset_temp)) - 2 # 2 == tri + freq_trinuc32 columns
  row_high = c(rep("high", ncols), "A(C>A)A")
  row_low = c(rep("low", ncols), "A(C>G)A")
  col_tri = data.frame(unique(offset_temp$tri)) %>% `colnames<-`("tri") %>% filter(!tri %in% c("A(C>A)A", "A(C>G)A"))
}


empty_offset = data.frame(matrix(ncol = length(names(offset_temp)))) %>% 
  `colnames<-`(names(offset_temp)) %>% 
  select(-c(freq_trinuc32)) %>% 
  rbind(row_high) %>% 
  rbind(row_low) %>% 
  drop_na %>% 
  {if (exists("col_repliseq")) merge(., col_repliseq, all = T) else .} %>% 
  merge(col_tri, all = T) %>% 
  expand(!!! syms(names(offset_temp)[!str_detect(names(offset_temp), "freq_trinuc32")])) %>%
  drop_na

if(!is.null(empty_offset$CTCF_cohesin)){
  empty_offset = empty_offset %>% 
    mutate(CTCF_cohesin = gsub("low", "CTCF_cohesin_peak", CTCF_cohesin),
           CTCF_cohesin = gsub("high", "bgGenome", CTCF_cohesin))
}

if(!is.null(empty_offset$A3A_TpCpH_hairpins)){
  empty_offset = empty_offset %>% 
    mutate(A3A_TpCpH_hairpins = gsub("low", "hairpin_TpCpH", A3A_TpCpH_hairpins),
           A3A_TpCpH_hairpins = gsub("high", "bgGenome", A3A_TpCpH_hairpins))
}

## correct feature levels (in case they are not "high" or "low")
# WARNING: this assumes that all features have 2 (and only 2) levels
for(feature in names(empty_offset)){
  empty_offset = empty_offset %>%
    mutate(!!sym(feature) := gsub("low", unique(offset_temp[[feature]])[1], !!sym(feature)),
           !!sym(feature) := gsub("high", unique(offset_temp[[feature]])[2], !!sym(feature)))
}
  
offset = merge(offset_temp, empty_offset, all = T) %>% 
  replace_na(list(freq_trinuc32 = 0)) %>% 
  mutate(trinuc32 = gsub("\\(", "", tri),
         trinuc32 = gsub("\\>.*\\)", "", trinuc32)) %>% 
  relocate(trinuc32, .before = tri)

write_tsv(offset, "offset.tsv")


###############################################
### QC
# check that distribution of trinucs after merging all features and chromosomes remain matched

# central base is bold
trinuc32_sorted = c("A𝗖A","A𝗖C","A𝗖G","A𝗖T","C𝗖A","C𝗖C","C𝗖G","C𝗖T","G𝗖A","G𝗖C","G𝗖G","G𝗖T","T𝗖A","T𝗖C","T𝗖G","T𝗖T","A𝗧A","A𝗧C","A𝗧G","A𝗧T","C𝗧A","C𝗧C","C𝗧G","C𝗧T","G𝗧A","G𝗧C","G𝗧G","G𝗧T","T𝗧A","T𝗧C","T𝗧G","T𝗧T")

trinuc_dist_QC = offset %>% 
  select(-tri) %>% 
  distinct %>% 
  pivot_longer(cols = !(matches("trinuc32") | matches("freq_trinuc32")),
               names_to = "feature", values_to = "levels") %>%
  group_by(trinuc32, feature, levels) %>%
  summarise(freq_trinuc32 = sum(freq_trinuc32)) %>%
  group_by(feature, levels) %>%
  mutate(prop_trinuc32 = freq_trinuc32 / sum(freq_trinuc32)) %>% 
  group_by(trinuc32, feature) %>%
  mutate(diff_prop_trinuc32 = abs(diff(prop_trinuc32)),
         diff_prop_trinuc32 = ifelse(diff_prop_trinuc32 * 10 >= 1,
                                            "!",
                                            ifelse(diff_prop_trinuc32 * 100 >= 1,
                                                   "***",
                                                   ifelse(diff_prop_trinuc32 * 1000 >= 1,
                                                          "**",
                                                          ifelse(diff_prop_trinuc32 * 10000 >= 1,
                                                                 "*",
                                                                 "")))),
         # keep the asterisks only in the lowest prop value of the 2 levels for each pair, to show in the plot cleanlier
         diff_prop_trinuc32 = ifelse(prop_trinuc32 == min(prop_trinuc32),
                                     diff_prop_trinuc32,
                                     "")) %>% 
  ungroup() %>%
  # central base is bold
  mutate(trinuc32 = gsub("C", "𝗖", trinuc32),
         trinuc32 = gsub("^𝗖", "C", trinuc32),
         trinuc32 = gsub("𝗖$", "C", trinuc32),
         trinuc32 = gsub("T", "𝗧", trinuc32),
         trinuc32 = gsub("^𝗧", "T", trinuc32),
         trinuc32 = gsub("𝗧$", "T", trinuc32),
         trinuc32 = factor(trinuc32, levels = trinuc32_sorted)) %>%
  mutate(levels = ifelse(levels %in% c("low", "0", "0-3", "0&1", "bgGenome"),
                         "low",
                         ifelse(levels %in% c("high", "1-3", "4-6", "2&3", "exon__non_CEG_TSG_OG"),
                                "high",
                                "OTHER LEVEL"))) %>%
  arrange(trinuc32, feature)

trinuc_dist_QC_plot = ggplot(trinuc_dist_QC,
                             aes(x = trinuc32,
                                 y = prop_trinuc32,
                                 fill = levels)) +
  scale_y_continuous(breaks = round(c(seq(0, 0.02, 0.01), 1/32, seq(0.04, max(trinuc_dist_QC$prop_trinuc32), 0.01)),2)) +
  geom_col(width = 0.6,
           position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(label = diff_prop_trinuc32),
            position = position_dodge(width = 0.75),
            hjust = -0.5, vjust = 0.8, size = 3, angle = 90) +
  geom_hline(yintercept = 1/32, linetype = "dashed") +
  facet_wrap(facets = ~feature,
             scales = "free") +
  xlab("32 trinucleotide types (collapsed from 64)") +
  ylab("Fraction of the 32 trinucleotide types") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 0),
        legend.position = c(0.85, 0),
        legend.justification = c(1, 0))

ggsave("trinuc_dist_QC.jpg",
       plot = trinuc_dist_QC_plot,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600)
