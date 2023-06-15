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
# Marcel K562 24 MMR-/- pairs
# Zou iPSC MMR-/-
# tumors MMRdef (hyper or MSI)

metadata_K562_MMRko_pairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv")
K562_MMRko_pairs = metadata_K562_MMRko_pairs %>% 
  select(`barcode sample`) %>% 
  `colnames<-`("sample_id")
write_tsv(K562_MMRko_pairs, "Marcel_K562_24_MMRko_pairs.tsv", col_names = F)

zou_MMRko = read_tsv("Zou_iPSC_MMRko.tsv", col_names = F) %>% 
  `colnames<-`("sample_id")

metadata_tumors_MMRdef = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv")
tumors_MMRdef = metadata_tumors_MMRdef %>% 
  filter(MSI_status %in% c("HYPER", "MSI")) %>% 
  select(sample_id)
write_tsv(tumors_MMRdef, "tumors_MMRdef.tsv")


### bind and write
write_csv(bind_rows(K562_MMRko_pairs,
                    zou_MMRko,
                    tumors_MMRdef), 
          "sample_ids.csv", col_names = F)
