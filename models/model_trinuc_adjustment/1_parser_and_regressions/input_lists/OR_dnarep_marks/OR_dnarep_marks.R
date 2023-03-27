library(tidyverse)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


dnarep_marks = read_csv("../dnarep_marks.csv")

offset = read_tsv(Sys.glob("../../work/[[:alnum:]][[:alnum:]]/*/offset.tsv")[1]) %>% 
  # revert log thing
  mutate(freq_trinuc32 = exp(log_freq_trinuc32) - 1) %>% 
  # keep cols of interest
  select(all_of(dnarep_marks$name), freq_trinuc32)


res_list = list()
dnarep_marks_left = dnarep_marks$name

for(ref_mark in dnarep_marks$name){
  
  if(length(dnarep_marks_left) <= 1){
    stop("Finished!\n")
    }
    
  dnarep_marks_left = str_remove(dnarep_marks_left, ref_mark)
  dnarep_marks_left = dnarep_marks_left[dnarep_marks_left != ""]
  
  or_table = offset %>% 
    select(all_of(ref_mark), all_of(dnarep_marks_left), freq_trinuc32) %>% 
    # "AID_target" to "high" since AID_target is alphabetically before bgGenome, so we keep the order high to low as in the remaining marks; also, makes sense, since AID_target means that the region is potentially "high" in AID mutagenesis
    mutate(!!ref_mark := gsub("AID_target", "high", get(ref_mark)),
           !!ref_mark := gsub("bgGenome", "low", get(ref_mark))) %>% 
    pivot_longer(cols = all_of(dnarep_marks_left), names_to = "grouping_mark", values_to = "level") %>% 
    group_by_at(vars(!matches("^freq_trinuc32$"))) %>% 
    summarise_at(vars(freq_trinuc32), ~sum(.)) %>% 
    pivot_wider(names_from = eval(ref_mark), values_from = freq_trinuc32) %>% 
    unite("grouping_mark_level", grouping_mark, level, sep = "__", remove = F) %>%
    select(-level) %>% relocate(grouping_mark) %>% 
    # arranging by 'grouping_mark_level' ensures that 'high' will be before 'low' for each grouping mark
    arrange(grouping_mark, grouping_mark_level)

  res_or = or_table %>% 
    # 'matrix(c(high, low)' ensures that first 'high' and then 'low', so the 1,1 --> 2,2 diagonal in each contingency table will be high×high and low×low
    summarise(OR = fisher.test(matrix(c(high, low), nrow = 2))$estimate)
  res_or_ci95_low = or_table %>% 
    summarise(CI95.low = fisher.test(matrix(c(high, low), nrow = 2))$conf.int[1])
  res_or_ci95_high = or_table %>% 
    summarise(CI95.high = fisher.test(matrix(c(high, low), nrow = 2))$conf.int[2])
  # res_fisher_pval = or_table %>% 
  #   summarise(fisher.pval = fisher.test(matrix(c(high, low), nrow = 2))$p)
  
  res = merge(res_or, res_or_ci95_low) %>% 
    merge(res_or_ci95_high) %>% #merge(res_fisher_pval)
    mutate(ref_mark = ref_mark)

  final_table = or_table %>% 
    mutate(grouping_mark_level = gsub(".*__", "", grouping_mark_level)) %>% 
    ungroup %>% 
    mutate(grouping_mark_level = gsub("AID_target", "high", grouping_mark_level),
           grouping_mark_level = gsub("bgGenome", "low", grouping_mark_level),
           grouping_mark_level = gsub("^", "grouping-", grouping_mark_level)) %>% 
    pivot_longer(cols = c(high, low), names_to = ref_mark, values_to = 'SNVs') %>% 
    mutate(!!ref_mark := gsub("^", "ref-", get(ref_mark))) %>% 
    unite("cell", grouping_mark_level, ref_mark) %>% 
    pivot_wider(names_from = cell, values_from = SNVs) %>% 
    left_join(res) %>% 
    relocate(ref_mark, .after = grouping_mark)
  
  res_list[[ref_mark]] = final_table
  gc()
}

res_table = bind_rows(res_list)

write_tsv(res_table, "OR_dnarep_marks.tsv")
#res_table = read_tsv("OR_dnarep_marks.tsv")


#######
## Heatmap

heatmap_table = res_table %>% 
  select(grouping_mark, ref_mark, OR) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% 
  pivot_wider(names_from = ref_mark, values_from = OR)

colname = names(heatmap_table)[2]
missing_row = heatmap_table[1:2] %>%
  pivot_wider(names_from = grouping_mark, values_from = colname) %>%
  mutate(!!colname := NA,
         grouping_mark = colname) %>%
  select(all_of(names(heatmap_table)))

heatmap_table = heatmap_table %>%
  bind_rows(missing_row)

missing_col = heatmap_table[1,] %>%
  select(-grouping_mark) %>%
  mutate(!!as.character(heatmap_table[1,1]) := NA) %>%
  pivot_longer(cols = everything(), names_to = "grouping_mark", values_to = as.character(heatmap_table[1,1]))

heatmap_table = heatmap_table %>%
  merge(missing_col, all = T) %>%
  mutate(!!colname := NA) %>%
  relocate(colname, .after = last_col()) %>%
  mutate(!!names(missing_col)[2] := ifelse(grouping_mark != colname,
                                           NA,
                                           get(names(missing_col)[2]))) %>%
  column_to_rownames("grouping_mark") %>% 
  arrange(desc(rowSums(is.na(.)))) %>% 
  as.matrix %>% 
  forceSymmetric(uplo="L") %>% 
  as.matrix


pdf(file = "OR_dnarep_marks.pdf", width = 10, height = 5.6)
heatmap = Heatmap(heatmap_table,
                  column_title = "      Odds-ratios of the genome-wide coordinates overlap between bins defined by their 'high' vs. 'low' abundance of a DNA-repair mark",
                  column_title_gp = gpar(fontsize = 10),
                  cell_fun = function(j, i, x, y, width, height, fill)
                  {
                    if (is.na(heatmap_table[i, j]) != TRUE){
                      grid.text(sprintf("%s", heatmap_table[i, j]), x, y, gp = gpar(fontsize = 10))
                    } else {
                      grid.text(sprintf(''), x, y, gp = gpar(fontsize = 10))
                    }
                  },
                  col = colorRamp2(c(min(heatmap_table, na.rm = T), 
                                     1,
                                     # 12 is too far away from next highest value, use 2nd highest value as red
                                     max(heatmap_table[heatmap_table != max(heatmap_table, na.rm = T)], na.rm = T)),
                                   c("blue", "white", "red")), 
                  na_col = "darkred",
                  border = "black",
                  show_column_names = T, 
                  show_row_names = T, 
                  cluster_rows = T, 
                  cluster_columns = T,
                  show_heatmap_legend = T,
                  heatmap_legend_param = list(title = "O.R.",
                                              legend_height = unit(10, "cm")),
                  row_names_rot = 0,
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize = 8),
                  column_names_rot = 45,
                  column_names_gp = gpar(fontsize = 8),
                  use_raster = T,
                  raster_by_magick = requireNamespace("magick", quietly = TRUE),
                  #raster_device = "tiff",
                  raster_quality = 600)
draw(heatmap)
dev.off()
