
### list with 100 samples that have more A3A T(C>D)H pairs closer than 1kb

# top 50 with morepairs --> A3A mutagenesis
top50_samples_morepairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv") %>% 
  arrange(desc(npairs)) %>% 
  slice_head(n = 50) %>% 
  arrange(desc(npairs)) %>% 
  select(sample_id)

# bottom with 0 (no A3A mutagenesis)
bottom50_samples_0pairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv") %>%  
  filter(npairs == 0) %>% 
  #slice_tail(n = 50) %>% 
  select(sample_id)

write_csv(bind_rows(top50_samples_morepairs,
                    bottom50_samples_0pairs), 
          "sample_ids.csv", col_names = F)
