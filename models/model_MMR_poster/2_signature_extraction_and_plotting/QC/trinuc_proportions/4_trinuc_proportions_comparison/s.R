offset_matched = read_tsv("../trinuc_proportions_matched_feature_tracks/offset.tsv")
offset_unmatched = read_tsv("../trinuc_proportions_unmatched_feature_tracks/res/offset.tsv")
# central base is bold
trinuc32_sorted = c("A𝗖A","A𝗖C","A𝗖G","A𝗖T","C𝗖A","C𝗖C","C𝗖G","C𝗖T","G𝗖A","G𝗖C","G𝗖G","G𝗖T","T𝗖A","T𝗖C","T𝗖G","T𝗖T","A𝗧A","A𝗧C","A𝗧G","A𝗧T","C𝗧A","C𝗧C","C𝗧G","C𝗧T","G𝗧A","G𝗧C","G𝗧G","G𝗧T","T𝗧A","T𝗧C","T𝗧G","T𝗧T")

trinuc_df_list = list()
for(trinuc in c("matched", "unmatched")){
  trinuc_df_list[[trinuc]] = get(paste0('offset_', trinuc)) %>% 
    select(-tri) %>% 
    distinct %>% 
    pivot_longer(cols = !(matches("trinuc32") | matches("freq_trinuc32")),
                 names_to = "feature", values_to = "levels") %>%
    group_by(trinuc32, feature, levels) %>%
    summarise(freq_trinuc32 = sum(freq_trinuc32)) %>%
    group_by(feature, levels) %>%
    mutate(prop_trinuc32 = freq_trinuc32 / sum(freq_trinuc32)) %>% 
    group_by(trinuc32, feature) %>%
    mutate(diff_prop_trinuc32 = abs(diff(prop_trinuc32))) %>% 
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
                                  "OTHER LEVEL")),
           matching = trinuc) %>%
    arrange(trinuc32, feature)
}

trinuc_dist_QC = bind_rows(trinuc_df_list) %>% 
  select(-c(levels, freq_trinuc32,prop_trinuc32)) %>% 
  distinct() %>% 
  arrange(trinuc32, feature)

# remove features that are not in the matched results (only in the unmatched)
features_to_rm = data.frame(table(trinuc_dist_QC$feature) - mean(table(trinuc_dist_QC$feature))) %>% filter(Freq<0) %>% pull(Var1) %>% as.character
trinuc_dist_QC = trinuc_dist_QC %>% 
  filter(!feature %in% features_to_rm) %>% 
  # x and y columns
  pivot_wider(names_from = matching, values_from = diff_prop_trinuc32)

trinuc_dist_QC_plot = ggplot(trinuc_dist_QC,
                             aes(x = unmatched,
                                 y = matched)) +
  geom_point() +
  geom_abline(slope = 1, intercept = c(0,0), linetype = "dashed", color = "darkgray") +
  ggrepel::geom_text_repel(aes(label = trinuc32),
                           max.overlaps = 10000,
                           size = 2,
                           force = 10) +
  facet_wrap(facets = ~feature,
             scales = "free") +
  ggtitle("Difference of each trinucleotide (32) proportions in 'high' vs. 'low' feature abundance") +
  xlab("Before matching") +
  ylab("After overlapping the coordinates of the matched features") +
  theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))

ggsave("trinuc_dist_comparison_QC.jpg",
       plot = trinuc_dist_QC_plot,
       device = "jpg",
       width = 20,
       height = 11,
       dpi = 600)
