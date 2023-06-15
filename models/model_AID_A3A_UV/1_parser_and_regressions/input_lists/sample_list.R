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



########################################

# 10 lymphoid with most SBS84exp
# 10 non-skin with most T(C>D)H-pairs
# 10 skin with most CC>TT
# 3 RPE1 with POLH-/- + UVC
# >80 PSCi (diverse treatments or gene-KOs)


### SBS84 exposure in lymph samples
top10_SBS84exp_lymph = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/AID_SHM/3_extract_SBS84__sort_AID_vs_nonAID_samples/lymphoid_samples_SBS84_exposure_K3.tsv") %>% 
  separate(Sample, into = c("rm1", "sample_id"), extra = "merge") %>% 
  mutate(sample_id = gsub(".*_", "", sample_id)) %>% 
  left_join(names_conversion_table) %>%
  mutate(sample_id2 = ifelse(is.na(sample_id2),
                              sample_id,
                              sample_id2)) %>% 
  select(sample_id2, SBS84_raw_exp_AIDtarg_sig1) %>% 
  rename("sample_id" = "sample_id2") %>% 
  # top10_SBS84
  arrange(desc(SBS84_raw_exp_AIDtarg_sig1)) %>% 
  slice_head(n = 10) %>% 
  select(sample_id)


### A3A pancancer (non-skin) 1kbpairs_TpC2DpH
top10_A3A_1kbpairs_TpC2DpH = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_count_1kbpairs_TpC2DpH_in_tumors/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv") %>% 
  # top10_A3A
  arrange(desc(TpC2DpH_pairs)) %>% 
  slice_head(n = 10) %>% 
  select(sample_id)


### SKIN samples (for CTCF UV)
samples_in_data_folder = read_csv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/sample_ids_tables/sample_ids_in_data_folder.csv", col_names = F)$X1
samples_skin = read_tsv("../../../../../metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv") %>% 
  filter(sample_id %in% samples_in_data_folder | sample_id_2 %in% samples_in_data_folder) %>% 
  # keep ONLY skin samples (to ensure that UV was acting on 5'-CC-3')
  filter(tissue == "skin" | str_detect(OriginalType, "[M,m]elanoma"))

top10_UV_skin = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/2_count_CC2TT_genomewide_in_tumors/res/samples_CC2TT_genomewide_counted.tsv") %>% 
  filter(sample_id %in% samples_skin) %>% 
  # top10_UV
  arrange(desc(CC2TT)) %>% 
  slice_head(n = 10) %>% 
  select(sample_id)


### POLH-/- UVC, iPSC Kucab, and iPSC Zou

top_cell_lines = read_tsv("cell_lines_samples.tsv")



### bind and write

write_csv(bind_rows(top10_SBS84exp_lymph,
                    top10_A3A_1kbpairs_TpC2DpH,
                    top10_UV_skin,
                    top_cell_lines), 
          "sample_ids.csv", col_names = F)
