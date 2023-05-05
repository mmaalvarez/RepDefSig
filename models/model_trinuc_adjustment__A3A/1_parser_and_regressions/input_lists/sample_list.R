
### calc ratio between 1kbpairs_TpC2DpH and total SNVs
A3A_1kbpairs_TpC2DpH_ratio = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv") %>% 
  mutate(ratio = A3A_1kbpairs_TpC2DpH/total_n_SNVs)


### list with 100 samples: 50 that have the most A3A mutagenesis and 50 that have the least

# top 50 with more A3A T(C>D)H pairs closer than 1kb (proportional to the total mut burden) --> A3A mutagenesis
top50_samples_morepairs = A3A_1kbpairs_TpC2DpH_ratio %>% 
  arrange(desc(ratio)) %>% 
  slice_head(n = 50) %>% 
  arrange(desc(ratio)) %>% 
  select(sample_id)

# top 50 with 0 1kbpairs_TpC2DpH AND the highest mut burden --> no A3A mutagenesis
bottom50_samples_0pairs = A3A_1kbpairs_TpC2DpH_ratio %>%  
  filter(A3A_1kbpairs_TpC2DpH == 0) %>% 
  arrange(desc(total_n_SNVs)) %>% 
  slice_head(n = 50) %>% 
  select(sample_id)

write_csv(bind_rows(top50_samples_morepairs,
                    bottom50_samples_0pairs), 
          "sample_ids.csv", col_names = F)
