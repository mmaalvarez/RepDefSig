library(tidyverse)
library(pROC)
library(data.table)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


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


##### dMMR vs MMRwt
K562 = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv") %>% 
  select(sample_id, `MMR deficiency expected`) %>% 
  rename("dMMR" = "MMR deficiency expected") %>% 
  mutate(group = "K562")
iPSC = read_tsv("../../../1_parser_and_regressions/input_lists/Zou_iPSC_MMRwt.tsv", col_names = F) %>% 
  mutate(dMMR = "MMRwt") %>% 
  bind_rows(read_tsv("../../../1_parser_and_regressions/input_lists/Zou_iPSC_MMRko.tsv", col_names = F) %>% 
              mutate(dMMR = "MMR-/-")) %>% 
  rename("sample_id" = "X1") %>% 
  mutate(group = "iPSC")
tumors = read_tsv("../../../1_parser_and_regressions/input_lists/tumors_MSI_MSS_metadata.tsv") %>% 
  replace_na(list(MSI_status = "NA")) %>% 
  rename("dMMR" = "MSI_status",
         "group" = "source") %>% 
  select(sample_id, dMMR, group)

samples_dMMR_status = bind_rows(K562, iPSC, tumors) %>% 
  rename("sample_id2" = "sample_id") %>% 
  left_join(names_conversion_table) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            sample_id2,
                            sample_id)) %>% 
  select(-sample_id2)


### load real-sample regression results
reg_coeffs_repdefsig = read_tsv("../../1_parser_and_regressions/res/results_real_samples.tsv") %>% #"../../../../../resources/OLD/OLD_MMR_poster/res/results_real_samples.tsv") %>% 
  rename_with(~str_replace(., 'high$', '')) %>% 
  rename_with(~str_replace(., 'bgGenome$', '')) %>% 
  rename("sample_id2" = "sample_id") %>% 
  left_join(names_conversion_table) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            sample_id2,
                            sample_id)) %>% 
  select(-c(info1,info2)) %>% 
  drop_na


###########################################################
### highlight (and then convert their -culprit- estimates and CIs to 0 for the PCA) samples with any estimate(s) >~|12|, as these would mean that the regression died
reg_coeffs_repdefsig = reg_coeffs_repdefsig %>% 
  mutate(betas_12 = ifelse(if_any(starts_with("estimate"), 
                                  ~abs(.x) >= 12), "yes", "no"),
         # convert the estimates (the ones that cause that betas_12 == yes) and their CI to 0, for PCA
         across(starts_with("estimate") | starts_with("conf."),
                ~ifelse(abs(.x) >= 12, 0, .)))
###########################################################


### PCA to look for biases regarding which samples were regressed with which theta values (either converged, or not)
pca_res = reg_coeffs_repdefsig %>% 
  select(sample_id, contains("estimate_")) %>% 
  column_to_rownames("sample_id") %>% 
  rename_with(~str_replace(., 'estimate_', '')) %>% 
  PCA(graph = FALSE)

pca = reg_coeffs_repdefsig %>%
  merge(pca_res$ind$coord %>% 
          data.frame %>% 
          rownames_to_column("sample_id")) %>% 
  select(sample_id, contains("Dim"), theta, betas_12) %>% 
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  mutate(sample_id = gsub("_REP1", "", sample_id)) %>% 
  left_join(samples_dMMR_status) %>% 
  rename("Sample" = "sample_id") %>% 
  unite("Source / dMMR status", group, dMMR, sep = " ")

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

outliers_cutoff = 1.4

pca_plot = ggplot(pca,
                  aes(x = PC1,
                      y = PC2)) +
  #coord_fixed() +
  geom_point(aes(size = theta,
                 fill = `Source / dMMR status`),
             shape = 21) +
  # highlight samples with any beta with coeff >|12| (reg. died)
  geom_point(data = filter(pca, betas_12 == "yes"),
             fill = "black",
             shape = 4,
             size = 5) +
  scale_fill_manual(values = jet.colors(length(unique(pca$`Source / dMMR status`)))) +
  guides(fill = guide_legend(override.aes = list(size=6, shape=21))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  ggrepel::geom_text_repel(data = pca[abs(pca$PC1) > outliers_cutoff*IQR(pca$PC1) |  abs(pca$PC2) > outliers_cutoff*IQR(pca$PC2),],
                           aes(label = Sample), # info1
                           size = 4,
                           force = 20,
                           segment.size = 0.1,
                           min.segment.length = 0.001,
                           max.overlaps = 100000,
                           max.iter = 100000) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent"),
        text = element_text(size = 25))
ggsave("pca.jpg",
       plot = pca_plot,
       device = "jpg",
       width = 20,
       height = 10,
       dpi = 600,
       bg = "white")

scree = fviz_eig(pca_res,
                 ylim = c(0, round(max(pca_res$eig[,2]))),
                 geom = c("bar"),
                 addlabels = FALSE,
                 ncp = length(rownames(pca_res$eig)),
                 main = "",
                 ggtheme = theme_classic(base_size = 20),
                 xlab = "PC",
                 ylab = "% variance")
ggsave("scree.jpg",
       plot = scree,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")

vars = pca_res$var$coord %>%
  data.frame() %>%
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  rownames_to_column("vars") %>% 
  ggplot() +
  geom_segment(aes(x = 0, xend = PC1,
                   y = 0, yend = PC2),
               arrow = arrow(length = unit(0.025,
                                           "npc"),
                             type = "open"),
               lwd = 0.5,
               linetype = "dashed") +
  ggrepel::geom_text_repel(aes(x = PC1,
                               y = PC2,
                               label = vars),
                           size = 6,
                           direction = "y",
                           vjust = 3,
                           force = 5,
                           segment.size = 0,
                           min.segment.length = 0,
                           max.iter = 100000) +
  theme_nothing()
ggsave("vars.jpg",
       plot = vars,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")



##########################
### only K562

pca_K562 = pca %>% 
  filter(str_detect(`Source / dMMR status`, "K562"))

pca_K562_plot = ggplot(pca_K562,
                  aes(x = PC1,
                      y = PC2)) +
  #coord_fixed() +
  geom_point(aes(size = theta,
                 fill = `Source / dMMR status`),
             shape = 21) +
  scale_fill_manual(values = jet.colors(length(unique(pca_K562$`Source / dMMR status`)))) +
  guides(fill = guide_legend(override.aes = list(size=6, shape=21))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  ggrepel::geom_text_repel(aes(label = Sample), # info1
                           size = 4,
                           force = 20,
                           segment.size = 0.1,
                           min.segment.length = 0.001,
                           max.overlaps = 100000,
                           max.iter = 100000) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent"),
        text = element_text(size = 25))
ggsave("pca_K562.jpg",
       plot = pca_K562_plot,
       device = "jpg",
       width = 20,
       height = 10,
       dpi = 600)



# ####################################################
# 
# #### AUROC and Mann whitney
# 
# 
# ### load SBS84-exp and regression results tables, and merge them
# 
# samples_with_high_exposure_other_than_SBS84 = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/AID_SHM/6_exploratory_plot/ground_truth_SBS84_exposures.tsv") %>% 
#   filter(group == "Lymphoid tumors") %>% 
#   pivot_longer(cols = starts_with("SNVs_"), names_to = "Genomic region", values_to = "SNVs") %>% 
#   mutate(`Genomic region` = gsub("SNVs_", "", `Genomic region`),
#          `Genomic region` = gsub("AID", "AID target ", `Genomic region`),
#          `Genomic region` = gsub("non", "non ", `Genomic region`),
#          `Genomic region` = factor(`Genomic region`, ordered = T, levels = c("non AID target regions", 
#                                                                              "AID target regions")),
#          ### WARNING -- ADD 1e-17 SNV to non-lymphoid that have 0 SNVs, so log10(SNVs) != Inf
#          SNVs = ifelse(group == "Non-lymphoid tumors" & SNVs == 0,
#                        1e-17,
#                        SNVs)) %>% 
#   pivot_longer(cols = contains("raw_exp"), names_to = "signature", values_to = "Raw exposure") %>% 
#   filter(!is.na(`Raw exposure`)) %>% 
#   mutate(signature = gsub("_raw_exp", "", signature)) %>% 
#   separate(signature, into = c("SBS-like", "Genomic_region_2", "signature"), sep = "_") %>% 
#   # WARNING SBS-like indicates the genome region that signature is based on (either AIDtargets or bgGenome)
#   unite("SBS-like", `SBS-like`,`Genomic_region_2`, sep = "-") %>% 
#   # WARNING summing up exposures to SBS9-like signatures (1 and 3) in lymphoid samples' bgGenome
#   group_by(source, sample_id, sample_id2, group,Tissue,`Genomic region`,SNVs, `SBS-like`) %>% 
#   summarise(`Raw exposure` = sum(`Raw exposure`)) %>% 
#   ungroup() %>% 
#   select(-c(sample_id2, Tissue)) %>% 
#   filter( (str_detect(`SBS-like`, paste(c("SBS7a", "SBS7b", "SBS9", "SBS2$"), collapse = "|")) & `Raw exposure` > 10)
#           |
#             (!(str_detect(`SBS-like`, "SBS84")) & `Raw exposure` > 100)
#   ) %>% 
#   select(sample_id) %>% 
#   distinct %>% 
#   pull(sample_id)
# 
# ground_truth_SBS84_exposures = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/AID_SHM/3_extract_SBS84__sort_AID_vs_nonAID_samples/lymphoid_samples_SBS84_exposure_K3.tsv") %>% 
#   separate(Sample, into = c("source", "sample_id"), sep = "_", extra = "merge") %>% 
#   separate(sample_id, into = c("cancer_type", "sample_id"), sep = "\\.", fill = "left") %>% 
#   separate(sample_id, into = c("country", "sample_id"), sep = "_", fill = "left") %>% 
#   dplyr::select(-c(Tissue, cancer_type, country)) %>% 
#   left_join(names_conversion_table) %>% 
#   mutate(sample_id2 = ifelse(is.na(sample_id2),
#                              sample_id,
#                              sample_id2)) %>% 
#   select(source,sample_id,sample_id2,SNVs_AIDregions,SNVs_nonAIDregions,SBS84_raw_exp_AIDtarg_sig1) %>%          
#   rename("SBS84 raw exposure" = "SBS84_raw_exp_AIDtarg_sig1") %>% 
#   ## keep only samples without hig exposures to SBS other than SBS84
#   filter(sample_id %in% samples_with_high_exposure_other_than_SBS84 | sample_id2 %in% samples_with_high_exposure_other_than_SBS84)
# 
# 
# SBS84exposures_regcoeffs = merge(ground_truth_SBS84_exposures, 
#                                  select(reg_coeffs_repdefsig, sample_id, contains("estimate_"), starts_with("conf"), theta)) %>% 
#   relocate(sample_id2, .after = sample_id) %>% 
#   arrange(desc(`SBS84 raw exposure`)) %>% 
#   rowwise() %>% 
#   # discretize Model prediction: if estimate's CI95% is above 0 there is AID SHM (since ref=AIDreg vs. alt=bgGenome)
#   mutate(`Model prediction` = ifelse(conf.high_AID_regions<0,
#                                      yes = "AID-SHM samples",
#                                      no = "non-AID-SHM samples"),
#          `Model prediction` = factor(`Model prediction`, ordered = T, levels = c("non-AID-SHM samples", "AID-SHM samples")))
# 
# write_tsv(SBS84exposures_regcoeffs, "SBS84exposures_regcoeffs.tsv")
# writexl::write_xlsx(SBS84exposures_regcoeffs, "SBS84exposures_regcoeffs.xlsx")
# 
# 
# #############################
# #### AUROC
# #############################
# 
# jpeg("AUROC_SBS84_raw_exposures.jpeg", width = 1700, height = 1500, quality = 100)
# par(mar = c(4, 4, 4, 4)+.1,
#     cex = 5.2)
# auroc = roc(data = SBS84exposures_regcoeffs,
#             predictor = `SBS84 raw exposure`,
#             response = `Model prediction`, levels = c("non-AID-SHM samples", "AID-SHM samples"),
#             direction="<", # SBS84 values in non-AID-SHM samples are lower than in AID-SHM samples
#             ci = T)
# plot(auroc,
#      main = "Ranked SBS84 raw exposures\nrecovering model's AID-SHM classification",
#      lwd = 10,
#      identity.lwd = 10)
# text(0.3, 0.4, paste("AUROC: ", 
#                      round(auroc$auc, 2),
#                      paste("\n(95% CI ", paste0(round(auroc$ci, 2)[1], "-", round(auroc$ci, 2)[3], ")"))), 
#      cex = 1.6)
# dev.off()
# 
# #############################
# #### mann-whitney
# #############################
# MW_plot_log10 = ggplot(SBS84exposures_regcoeffs %>% 
#                          mutate(`Model prediction` = as.character(`Model prediction`),
#                                 `Model prediction` = ifelse(`Model prediction` == "non-AID-SHM samples",
#                                                             "non-AID-SHM samples\n(samples with y=0 not shown)", 
#                                                             `Model prediction`),
#                                 `Model prediction` = factor(`Model prediction`, ordered = T, 
#                                                             levels = c("non-AID-SHM samples\n(samples with y=0 not shown)", 
#                                                                        "AID-SHM samples"))),
#                        aes(x = `Model prediction`,
#                            y = `SBS84 raw exposure`)) +
#   scale_y_log10() +
#   coord_cartesian(ylim = c(1, max(SBS84exposures_regcoeffs$`SBS84 raw exposure`))) +
#   geom_violin(fill = NA,
#               linewidth = 2) +
#   geom_boxplot(fill = NA, linewidth = 2, width = 0.5, notch = T, outlier.alpha = 0) +
#   ggbeeswarm::geom_beeswarm(aes(fill = source),
#                             shape = 21,
#                             size = 4,
#                             alpha = 0.8) +
#   scale_fill_manual(values = c("red", "blue", "yellow")) +
#   guides(fill = guide_legend(override.aes = list(size = 5))) +
#   annotate("text", x = 1.5, y = 1000, size = 8,
#            label = paste0("Mann-Whitney 'less-than'\np-value = ",
#                           as.character(signif(wilcox.test(`SBS84 raw exposure` ~ `Model prediction`,
#                                                           SBS84exposures_regcoeffs,
#                                                           alternative = "less")[["p.value"]], digits = 2)))) + 
#   theme_bw() +
#   ggtitle("Lymphoid samples") +
#   xlab("Model prediction\n#SNVs ~ 9×DNArep + WG-vs-AIDtargets + RepliSeq + offset") +
#   ylab("SBS84 raw exposure (log10 scale)") +
#   theme(text = element_text(size = 25),
#         legend.title = element_blank())
# ggsave("MW_plot_SBS84_raw_exposures_log10.jpg",
#        plot = MW_plot_log10,
#        device = "jpg",
#        width = 19,
#        height = 10,
#        dpi = 600,
#        bg = "white")
