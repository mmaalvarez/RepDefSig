PCAWG_names_conversion_table = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/pre_metadata/OLD_metadata.tsv") %>%
  rename("sample_id2" = "sample_id",
         "sample_id" = "icgc_donor_id") %>%
  filter(source %in% c("PCAWG")) %>% 
  mutate(icgc_donor_id = NA) %>% 
  dplyr::select(icgc_donor_id, sample_id, sample_id2)

TCGA_names_conversion_table = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/sample_ids_tables/TCGA_sampleid_conversion.tsv")
full_TCGA_names = data.frame(icgc_donor_id = as.character(),
                             sample_id2 = as.character())
for(TCGA_sample in TCGA_names_conversion_table$sample_id){
  
  # find full sample name in ../../data/
  full_sample_name = list.files('/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data', pattern = TCGA_sample)
  
  if(length(full_sample_name) >= 1){
    new_entry = data.frame(icgc_donor_id = as.character(filter(TCGA_names_conversion_table, sample_id == TCGA_sample) %>% pull(icgc_donor_id)),
                           sample_id2 = as.character(full_sample_name)) %>% 
      mutate(sample_id2 = gsub("muts_pass_", "", sample_id2),
             sample_id2 = gsub(".csv", "", sample_id2))
    # append to table
    full_TCGA_names = bind_rows(full_TCGA_names, new_entry)
  }
}
TCGA_names_conversion_table = merge(TCGA_names_conversion_table, full_TCGA_names) %>% 
  dplyr::select(icgc_donor_id, sample_id, sample_id2)

names_conversion_table = bind_rows(PCAWG_names_conversion_table, TCGA_names_conversion_table) %>% 
  mutate(sample_id = ifelse(!is.na(icgc_donor_id),
                            icgc_donor_id,
                            sample_id)) %>% 
  dplyr::select(sample_id, sample_id2)



### SBS84 exposure in lymph samples
SBS84exp_lymph = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/AID_SHM/3_extract_SBS84__sort_AID_vs_nonAID_samples/lymphoid_samples_SBS84_exposure_K3.tsv") %>% 
  arrange(desc(SBS84_raw_exp_AIDtarg_sig1)) %>% 
  separate(Sample, into = c("rm1", "sample_id"), extra = "merge") %>% 
  mutate(sample_id = gsub(".*_", "", sample_id)) %>% 
  left_join(names_conversion_table) %>%
  mutate(sample_id2 = ifelse(is.na(sample_id2),
                              sample_id,
                              sample_id2)) %>% 
  select(sample_id2, SBS84_raw_exp_AIDtarg_sig1) %>% 
  rename("sample_id" = "sample_id2")

### A3A pancancer: calc ratio between 1kbpairs_TpC2DpH and total SNVs
A3A_1kbpairs_TpC2DpH_ratio = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv") %>% 
  mutate(ratio = A3A_1kbpairs_TpC2DpH/total_n_SNVs)


### create list with:
#             - 25 lymphoid that have the most SBS84 exposure (AID mutagenesis)
#             - 25 lymphoid that have the least SBS84 exposure
#             - 25 pancancer that have the most A3A mutagenesis: more A3A T(C>D)H pairs closer than 1kb (proportional to the total mut burden)
#             - 25 pancancer that have the least A3A mutagenesis: 0 1kbpairs_TpC2DpH AND the highest mut burden

top25_SBS84 = SBS84exp_lymph %>% 
  slice_head(n = 25) %>% 
  select(sample_id)
bottom25_SBS84 = SBS84exp_lymph %>% 
  slice_tail(n = 25) %>% 
  select(sample_id)
top25_A3A = A3A_1kbpairs_TpC2DpH_ratio %>% 
  arrange(desc(ratio)) %>% 
  slice_head(n = 25) %>% 
  select(sample_id)
bottom25_A3A = A3A_1kbpairs_TpC2DpH_ratio %>%  
  filter(A3A_1kbpairs_TpC2DpH == 0) %>% 
  arrange(desc(total_n_SNVs)) %>% 
  slice_head(n = 25) %>% 
  select(sample_id)

write_csv(bind_rows(top25_SBS84,
                    bottom25_SBS84,
                    top25_A3A,
                    bottom25_A3A), 
          "sample_ids.csv", col_names = F)
