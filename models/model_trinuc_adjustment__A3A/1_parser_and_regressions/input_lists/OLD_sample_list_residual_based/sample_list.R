
### list with 100 samples that have positive residuals:

# top 50 largest residuals --> A3A mutagenesis
top50_samples_pos_residuals = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/3_exploratory_plot_and_sample_list/residuals.tsv") %>% 
  filter(residuals >=-1) %>% 
  arrange(desc(residuals)) %>% 
  slice_head(n = 50) %>% 
  arrange(desc(residuals)) %>% 
  select(sample_id)

# bottom 50 smallest positive residuals (i.e. closer to 0 --> no A3A mutagenesis), allowing for some negative decimals below 0 (>=-0.03)
bottom50_samples_pos_residuals = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/3_exploratory_plot_and_sample_list/residuals.tsv") %>% 
  filter(residuals >=-3e-2) %>% 
  arrange(desc(residuals)) %>% 
  slice_tail(n = 50) %>% 
  arrange(desc(residuals)) %>% 
  select(sample_id)

write_csv(bind_rows(top50_samples_pos_residuals,
                    bottom50_samples_pos_residuals), 
          "sample_ids.csv", col_names = F)
