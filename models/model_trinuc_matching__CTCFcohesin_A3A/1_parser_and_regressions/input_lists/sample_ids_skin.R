samples_in_data_folder = read_csv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/sample_ids_tables/sample_ids_in_data_folder.csv", col_names = F)$X1

sample_ids_skin = read_tsv("../../../../../metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv") %>% 
  filter(sample_id %in% samples_in_data_folder | sample_id_2 %in% samples_in_data_folder) %>% 
  # keep ONLY skin samples (to ensure that UV was acting on 5'-CC-3')
  filter(tissue == "skin" | str_detect(OriginalType, "[M,m]elanoma")) %>% 
  select(sample_id)

write_tsv(sample_ids_skin,
          "sample_ids_skin.csv", 
          col_names = F)
