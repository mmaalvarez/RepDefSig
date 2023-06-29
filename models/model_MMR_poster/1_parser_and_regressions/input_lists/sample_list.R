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
# Marcel K562 24 KO pairs (strong, mid, weak, or no dMMR)
# Zou iPSC MMR-/- + controls
# tumors MMRdef (hyper or MSI) + MMRwt (MSS)


metadata_K562_MMRko_pairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv")
K562_MMRko_pairs = metadata_K562_MMRko_pairs %>% 
  select(`barcode sample`) %>% 
  `colnames<-`("sample_id")
write_tsv(K562_MMRko_pairs, "Marcel_K562_24_MMRko_pairs.tsv", col_names = F)



zou_MMRko = read_tsv("Zou_iPSC_MMRko.tsv", col_names = F) %>% 
  `colnames<-`("sample_id")
# 8 controls + 1 "MMRko" sample with very low mut burden (123_s1)
zou_MMRwt = read_tsv("Zou_iPSC_MMRwt.tsv", col_names = F) %>% 
  `colnames<-`("sample_id")



metadata_tumors = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv")

tumors_MSI = metadata_tumors %>% 
  filter(MSI_status %in% c("HYPER", "MSI"))

tumors_MSS = metadata_tumors %>% 
  ## keep things as similar as possible to the MSI samples (except the MSI status)
  filter(source %in% c("HMF", "TCGA", "PCAWG")) %>% 
  filter(tissue %in% unique(tumors_MSI$tissue)) %>% 
  # (before we were missing 1 HMF ovary, so specify that)
  filter(OriginalType %in% unique(tumors_MSI$OriginalType) | str_detect(OriginalType, "Ovary")) %>% 
  filter(hr_status %in% unique(tumors_MSI$hr_status)) %>%
  filter(smoking_history %in% unique(tumors_MSI$smoking_history)) %>%
  filter(treatment_platinum %in% unique(tumors_MSI$treatment_platinum)) %>%
  filter(treatment_5FU %in% unique(tumors_MSI$treatment_5FU)) %>%
  filter(gender %in% unique(tumors_MSI$gender)) %>% 
  ## specify only MSS (or NA in the case of TCGA, as it does not say MSS explicitly, it just says MSI or NA)
  filter(!MSI_status %in% c("HYPER", "MSI", "ERROR", "ERCC2mut")) %>% 
  filter(!is.na(MSI_status) | source == "TCGA")

## now randomly keep 162 MSS samples that keep tumors_MSI metadata proportions

# proportions of source and tissue to keep in MSS samples
tumors_MSI_metadata_proportions = tumors_MSI %>% 
  select(source, tissue) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq > 0) %>% 
  arrange(desc(Freq), source, tissue)

# function to apply sample_n for each group
sample_group <- function(group_df, freq) {
  group_df %>% 
    sample_n(min(freq, n()))
}

# apply
tumors_MSS = tumors_MSI_metadata_proportions %>% 
  pmap_dfr(~tumors_MSS %>%
             filter(source == ..1, tissue == ..2) %>%
             sample_group(freq = ..3))

# check that retained proportions are the same as for MSI
tumors_MSS_metadata_proportions = tumors_MSS %>% 
  select(source, tissue) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq > 0) %>% 
  arrange(desc(Freq), source, tissue)


# bind MSI and MSS
write_tsv(bind_rows(tumors_MSI, tumors_MSS), "tumors_MSI_MSS_metadata.tsv")



### bind all and write out
write_csv(bind_rows(K562_MMRko_pairs,
                    zou_MMRko,
                    zou_MMRwt,
                    select(tumors_MSI, sample_id),
                    select(tumors_MSS, sample_id)), 
          "sample_ids.csv", col_names = F)
