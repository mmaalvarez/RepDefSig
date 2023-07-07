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
K562_MMRko_pairs = metadata_K562_MMRko_pairs
write_tsv(K562_MMRko_pairs, "Marcel_K562_24_MMRko_pairs.tsv")



zou_MMRko = read_tsv("Zou_iPSC_MMRko.tsv", col_names = F) %>% 
  `colnames<-`("sample_id")
# 8 controls + 1 "MMRko" sample with very low mut burden (123_s1)
zou_MMRwt = read_tsv("Zou_iPSC_MMRwt.tsv", col_names = F) %>% 
  `colnames<-`("sample_id")



metadata_tumors = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv")

tumors_MSI = metadata_tumors %>% 
  filter((MSI_status %in% c("HYPER", "MSI") & (is.na(msStatus) | msStatus != "MSS")) |
         (msStatus == "MSI" & (is.na(MSI_status) | !MSI_status %in% c("MSS", "ERROR", "ERCC2mut"))))

tumors_MSS = metadata_tumors %>% 
  ## keep things as similar as possible to the MSI samples (except the MSI status)
  filter(source %in% unique(tumors_MSI$source)) %>% # it's "HMF", "TCGA", and "PCAWG" 
  filter(tissue %in% unique(tumors_MSI$tissue)) %>%
  ## specify only MSS (or NA in the case of TCGA, as it does not say MSS explicitly, it just says MSI or NA)
  filter((!MSI_status %in% c("HYPER", "MSI", "ERROR", "ERCC2mut") | is.na(MSI_status)) &
         (msStatus != "MSI" | is.na(msStatus))) %>% 
  filter(!(is.na(MSI_status) & is.na(msStatus)))

## now keep length(tumors_MSI$sample_id) MSS samples that keep tumors_MSI metadata proportions

# proportions of source and tissue to keep in MSS samples
tumors_MSI_metadata_proportions = tumors_MSI %>% 
  select(source, tissue) %>% 
  table %>% 
  data.frame %>% 
  filter(Freq > 0) %>% 
  arrange(desc(Freq), source, tissue)

# function to slice top n=='freq' samples with lowest msIndelsPerMp (i.e. least MSI) for each source-tissue group
sample_group <- function(group_df, freq) {
  group_df %>% 
    arrange(msIndelsPerMb) %>% 
    slice_head(n = freq)
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
write_csv(bind_rows(select(K562_MMRko_pairs, sample_id),
                    zou_MMRko,
                    zou_MMRwt,
                    ## tumors have been downsampled, so they have the '.downsampled_SNVs' suffix now
                    select(tumors_MSI, sample_id) %>% 
                      mutate(sample_id = gsub("$", ".downsampled_SNVs", sample_id)),
                    select(tumors_MSS, sample_id) %>% 
                      mutate(sample_id = gsub("$", ".downsampled_SNVs", sample_id))), 
          "sample_ids.csv", col_names = F)
