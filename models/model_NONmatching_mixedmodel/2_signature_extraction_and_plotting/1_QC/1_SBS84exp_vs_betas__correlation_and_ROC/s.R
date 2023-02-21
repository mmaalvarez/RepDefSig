library(tidyverse)
library(pROC)


## name conversion tables 

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



### load SBS84-exp and regression results tables, and merge them

ground_truth_SBS84_exposures = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/AID_SHM/3_extract_SBS84__sort_AID_vs_nonAID_samples/lymphoid_samples_SBS84_exposure_K3.tsv") %>% 
  separate(Sample, into = c("source", "sample_id"), sep = "_", extra = "merge") %>% 
  separate(sample_id, into = c("cancer_type", "sample_id"), sep = "\\.", fill = "left") %>% 
  separate(sample_id, into = c("country", "sample_id"), sep = "_", fill = "left") %>% 
  dplyr::select(-c(Tissue, cancer_type, country)) %>% 
  left_join(names_conversion_table) %>% 
  mutate(sample_id2 = ifelse(is.na(sample_id2),
                             sample_id,
                             sample_id2)) %>% 
  rename("SBS84 raw exposure" = "SBS84_raw_exp",
         "SBS84 normalized exposure" = "SBS84_norm_exp")

reg_coeffs_repdefsig = read_tsv("../1_parser_and_regressions/res/results.tsv") %>% 
  dplyr::select(-theta) %>% 
  rename_with(~str_replace(., '.L', '')) %>% 
  rename_with(~str_replace(., '_AID_regions', '')) %>% 
  rename("sample_id2" = "sample_id") %>% 
  left_join(names_conversion_table) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            sample_id2,
                            sample_id)) 

SBS84exposures_regcoeffs = left_join(ground_truth_SBS84_exposures, reg_coeffs_repdefsig) %>% 
  relocate(sample_id2, .after = sample_id) %>% 
  arrange(desc(`SBS84 raw exposure`)) %>% 
  rowwise() %>% 
  # discretize Model prediction: if estimate's CI95% is above 0 (and estimate is larger than estimates' median), there is AID SHM
  mutate(`Model prediction` = ifelse(conf.low>0 & conf.high>0 & estimate >= median(reg_coeffs_repdefsig$estimate),
                                     yes = "AID-SHM samples",
                                     no = "non-AID-SHM samples"),
         `Model prediction` = factor(`Model prediction`, ordered = T, levels = c("non-AID-SHM samples", "AID-SHM samples")))

write_tsv(SBS84exposures_regcoeffs, "SBS84exposures_regcoeffs.tsv")
writexl::write_xlsx(SBS84exposures_regcoeffs, "SBS84exposures_regcoeffs.xlsx")



#############################################
#### AUROC
jpeg("AUROC_SBS84_raw_exposures.jpeg", width = 1700, height = 1500, quality = 100)
par(mar = c(4, 4, 4, 4)+.1,
    cex = 5.2)
auroc = roc(data = SBS84exposures_regcoeffs,
            predictor = `SBS84 raw exposure`,
            response = `Model prediction`, levels = c("non-AID-SHM samples", "AID-SHM samples"),
            direction="<", # SBS84 values in non-AID-SHM samples are lower than in AID-SHM samples
            ci = T)
plot(auroc,
     main = "Ranked SBS84 raw exposures\nrecovering model's AID-SHM classification",
     lwd = 10,
     identity.lwd = 10)
text(0.3, 0.4, paste("AUROC: ", 
                     round(auroc$auc, 2),
                     paste("\n(95% CI ", paste0(round(auroc$ci, 2)[1], "-", round(auroc$ci, 2)[3], ")"))), 
     cex = 1.6)
dev.off()

#############################
#### mann-whitney
MW_plot_log10_no_2outliers = ggplot(SBS84exposures_regcoeffs %>% 
                                      mutate(`Model prediction` = as.character(`Model prediction`),
                                             `Model prediction` = ifelse(`Model prediction` == "non-AID-SHM samples",
                                                                         "non-AID-SHM samples\n(two samples with y=0 not shown)", 
                                                                         `Model prediction`),
                                             `Model prediction` = factor(`Model prediction`, ordered = T, 
                                                                         levels = c("non-AID-SHM samples\n(two samples with y=0 not shown)", 
                                                                                    "AID-SHM samples"))),
                                    aes(x = `Model prediction`,
                                        y = `SBS84 raw exposure`)) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1, max(SBS84exposures_regcoeffs$`SBS84 raw exposure`))) +
  geom_violin(fill = NA,
              linewidth = 2) +
  geom_boxplot(fill = NA, linewidth = 2, width = 0.5, notch = T, outlier.alpha = 0) +
  ggbeeswarm::geom_beeswarm(aes(fill = source),
                            shape = 21,
                            size = 4,
                            alpha = 0.8) +
  scale_fill_manual(values = c("red", "blue", "yellow")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  annotate("text", x = 1.5, y = 1000, size = 8,
           label = paste0("Mann-Whitney 'less-than'\np-value = ",
                          as.character(signif(wilcox.test(`SBS84 raw exposure` ~ `Model prediction`,
                                                          SBS84exposures_regcoeffs,
                                                          alternative = "less")[["p.value"]], digits = 2)))) + 
  theme_bw() +
  ggtitle("Lymphoid samples") +
  xlab("Model prediction\n#SNVs ~ WG-vs-AIDtargets + (1|RepliSeq) + (1|3NContext) + log(sum(RepliSeq×3NContext))") +
  ylab("SBS84 raw exposure (log10 scale)") +
  theme(text = element_text(size = 25),
        legend.title = element_blank())
ggsave("MW_plot_SBS84_raw_exposures_log10_no2outliers.jpg",
       plot = MW_plot_log10_no_2outliers,
       device = "jpg",
       width = 19,
       height = 10,
       dpi = 600,
       bg = "white")


##################################
### Exploratory plot

## NEW add controls (i.e. non-lymphoid tumor SBS84 exposures)

nonlymphoid = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/AID_SHM/5_non_lymphoid_controls/3_extract_SBS84__sort_AID_vs_nonAID_samples/nonlymphoid_samples_SBS84_exposure_K10.tsv")

SBS84exposures_regcoeffs = SBS84exposures_regcoeffs %>% 
  # fix this
  bind_rows(nonlymphoid)

exploratory_plot = ggplot(SBS84exposures_regcoeffs %>% 
         pivot_longer(cols = starts_with("SNVs_"), names_to = "Genomic region", values_to = "SNVs") %>% 
         mutate(`Genomic region` = gsub("SNVs_", "", `Genomic region`),
                `Genomic region` = gsub("AID", "AID target ", `Genomic region`),
                `Genomic region` = gsub("non", "non ", `Genomic region`),
                `Genomic region` = factor(`Genomic region`, ordered = T, levels = c("non AID target regions", 
                                                                                    "AID target regions"))),
       aes(x = `SBS84 raw exposure`,
           y = `SNVs`)) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(1, max(SBS84exposures_regcoeffs$`SBS84 raw exposure`)),
                  ylim = c(1, max(SBS84exposures_regcoeffs$SNVs_nonAIDregions))) +
  geom_line(aes(group = `SBS84 raw exposure`),
            alpha = 0.5,
            col = "darkgray",
            linetype = "dashed") +
  geom_point(aes(color = `Genomic region`,
                 shape = `Model prediction`),
             size = 3,
             alpha = 0.8) +
  # quadratic
  stat_smooth(aes(group = `Genomic region`),
              method = "lm",
              formula = y ~ poly(x, degree=2), size=1, col="black") +
  ggpmisc::stat_poly_eq(parse=T,
                        aes(group = `Genomic region`,
                            col = `Genomic region`,
                            label = ..eq.label..),
                        size = 7,
                        formula = y ~ poly(x, degree=2)) +
  scale_color_manual(values = c("blue", "red")) +
  scale_shape_manual(values = c(4, 19)) +
  xlab("SBS84 raw exposure\n(log10 scale; two samples with x=0 not shown)") +
  ylab("SNVs (log10 scale)") +
  theme_bw() +
  ggtitle("Lymphoid samples") +
  theme(#legend.title = element_blank(),
        #legend.position = "top",
        text = element_text(size = 25))
ggsave("exploratory_plot.jpg",
       plot = exploratory_plot,
       device = "jpg",
       width = 19,
       height = 10,
       dpi = 600,
       bg = "white")
