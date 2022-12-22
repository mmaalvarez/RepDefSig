library(tidyverse)
library(data.table)
library(msm) # rtnorm
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")



##### load and parse inputs

metadata = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x) %>% 
           select(sample_id, starts_with("info"))) %>% 
  rename("Sample" = "sample_id") %>% 
  mutate(dataset = ifelse(str_detect(Sample, "MSM0"),
                          "Kucab et al. 2019",
                          ifelse(str_detect(Sample, "MSK0"),
                                 "Zou et al. 2021",
                                 "ERROR: Unexpected sample name")),
         info1 = ifelse(str_detect(Sample, "MSK0"),
                        paste0(info1, "ko"),
                        ifelse(str_detect(Sample, "MSM0"),
                               info1,
                               "ERROR: Unexpected sample name")),
         info2 = gsub("^[a-z]_", "", info2), 
         info2 = gsub("DNA damage response inhibitors", "DNA damage resp. inh.", info2),
         info2 = ifelse(str_detect(Sample, "MSK0"),
                        paste0(info2, " (pathway)"),
                        ifelse(str_detect(Sample, "MSM0"),
                               paste0(info2, " (treatment)"),
                               "ERROR: Unexpected sample name")))

results_regressions = read_tsv("../../1_parser_and_regressions/res/results.tsv") %>% #../../1_parser_and_regressions/bin/output/results.tsv") %>% #../../../model1/1_parser_and_regressions/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_"), contains("conf"), glm)



#### Generate matrices resampling coefficients from their CI95% distributions (instead of UPmultinomial)

coefficient_table = results_regressions %>%
  select(sample_id, contains("estimate_"), contains("conf")) %>% 
  rename_all(~str_replace_all(., '.L', '')) %>% 
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat) %>% 
  group_by(sample_id, mark)

resample_from_CI = function(coefficient_table){
    summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                             mean = estimate,
                                                             sd = 1, 
                                                             lower = conf.low,
                                                             upper = conf.high))}
totalNumIters = 10000
coefficient_matrix_Resamp = list()

for (nIter in 1:totalNumIters) {
  cat(sprintf("Generating bootstrap matrix: nIter %d\n", nIter))
  
  # for each sample (row) resample coefficients from CI95% distrs.
  coefficient_matrix_TempIter = resample_from_CI(coefficient_table) %>% 
    pivot_wider(names_from = mark, values_from = resampled_estimate) %>% 
    ungroup %>% 
    column_to_rownames("sample_id") %>%
    data.matrix
  
  ## scale the CI-permuted coefficients matrices to [-1, +1] (e.g. python minmax scaler), because final decoding layer uses tanh as activation
  scale...

  coefficient_matrix_Resamp[[nIter]] = coefficient_matrix_TempIter
  gc()
}



################################
###### run the autoencoder



#################################



######
### prepare plot inputs

# Signature exposures in samples
exposures = autoencoder_res$h %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("Signature") %>%
  pivot_longer(cols = !contains("Signature"), names_to = "Sample", values_to = "Exposure") %>% 
  # add metadata info (e.g. treatments, MSI, HR, smoking...)
  left_join(metadata) %>% 
  mutate(Signature = factor(Signature, levels = unique(rownames(data.frame(autoencoder_res$h))))) %>%
  # highlight top hits for each signature
  group_by(Signature) %>% 
  mutate(is.hit = ifelse(Exposure==max(Exposure), "Top hit", NA))

# DNA repair mark weights in signatures
weights = autoencoder_res$w %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column("dna_repair_mark") %>%
  pivot_longer(cols = contains("autoencoder"), names_to = "signature", values_to = "weight") %>% 
  extract(dna_repair_mark, into = c("dna_repair_mark", "submatrix"), "(.*)_([^_]+$)") %>% 
  mutate(dna_repair_mark = sub(".L", "", dna_repair_mark),
         dna_repair_mark = factor(dna_repair_mark, levels = unique(gsub(".L_...coeff", "", rownames(autoencoder_res$w)))),
         signature = factor(signature, levels = unique(exposures$Signature))) %>% 
  arrange(dna_repair_mark, signature) %>% 
  group_by(dna_repair_mark, signature) %>% 
  ## subtract results of pos - results of neg, and get the absolute
  summarise(Weight = abs(.Primitive("-")(weight[1], weight[2]))) %>% 
  ungroup %>% 
  rename("DNA repair\nactivity" = "dna_repair_mark", 
         "Signature" = "signature") %>% 
  relocate(Signature)



########
## calculate "good-model-score": 1 would mean that ALL samples have 100% exposure from the signature(s) that is contributed 100% by a specific DNArep mark whose mechanism/pathway overlaps completely with the sample´s condition (gene-/-, treatment...), and 0 the opposite
# actually only considering "ground-truth" sample_condition--DNArep_pathway pairs ('condition_pathway_pairs')

repair_mark_pathways = results_regressions %>% 
  select(contains("estimate")) %>% 
  rename_with(~str_replace(., 'estimate_', '')) %>% 
  rename_with(~str_replace(., '.L', '')) %>% 
  colnames %>% 
  data.frame %>% 
  `colnames<-`("dna_repair_mark") %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(dna_repair_mark, "OGG1_"),
                                                  "BER",
                                                  ifelse(str_detect(dna_repair_mark, "UV_"),
                                                         "NER",
                                                         ifelse(str_detect(dna_repair_mark, "MSH6_|SETD2_"),
                                                                "MMR",
                                                                ifelse(str_detect(dna_repair_mark, "XRCC4"),
                                                                       "DSBR",
                                                                       NA)))))
condition_pathway_pairs = metadata %>% 
  select(info1, info2) %>% 
  distinct %>% 
  filter(! str_detect(info2, "[C,c]ontrol")) %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(info1, "[G,g]amma"),
                                                   "BER,DSBR",
                                                   ifelse(str_detect(info2, "BER |[A,a]lkylating|[N,n]itrosamine"),
                                                          "BER",
                                                          ifelse(str_detect(info2, "NER |[A,a]romatic|[H,h]eterocyclic|PAH|Radiation") |
                                                                   str_detect(info1, "platin"),
                                                                 "NER",
                                                                 ifelse(str_detect(info2, "MMR "),
                                                                        "MMR",
                                                                        ifelse(str_detect(info2, "DSB|DNA damage resp. inh.|[H,h]elicas|NHEJ|MMEJ|HR ") |
                                                                                 str_detect(info1, "PARP|Etoposide|Bleomycin|Camptothecin|Olaparib|Temozolomide|Melphalan|Cyclophosphamide|Mechlorethamine"),
                                                                               "DSBR",
                                                                               NA)))))) %>% 
  separate_rows(putative_repair_pathway_involved, sep = ",") %>% 
  filter(!is.na(putative_repair_pathway_involved)) %>% 
  left_join(metadata) %>% 
  filter(Sample %in% results_regressions$sample_id) %>% 
  arrange(putative_repair_pathway_involved) %>% 
  left_join(repair_mark_pathways) %>% 
  select(Sample, dna_repair_mark) %>% 
  rename("DNA repair\nactivity" = "dna_repair_mark") %>% 
  group_by(Sample) %>%
  mutate(Sample_possible_marks = n()) %>%
  ungroup()
condition_pathway_pairs = condition_pathway_pairs %>% 
  rowwise %>% 
  mutate(max_sample_mark_score = 1 / length(unique(condition_pathway_pairs$Sample)) / Sample_possible_marks) %>%
  select(-c(Sample_possible_marks))

good_model_score_table = weights %>% 
  merge(condition_pathway_pairs) %>% 
  merge(exposures) %>% 
  rename("sensical sample-mark pair" = "DNA repair\nactivity",
         "mark weight" = "Weight",
         "signature exposure" = "Exposure") %>% 
  select(Sample, "sensical sample-mark pair", Signature, "mark weight", "signature exposure", max_sample_mark_score) %>% 
  mutate(sample_mark_score = `mark weight` * `signature exposure` * max_sample_mark_score) %>% 
  arrange(Sample, `sensical sample-mark pair`, Signature)

write_tsv(good_model_score_table, "good_model_score_table.tsv")

good_model_score = sum(good_model_score_table$sample_mark_score)



###############
##### plotting

exposure_limit = 0.025

jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))

# plot only samples with Exposures > exposure_limit%
filtered_exposures = filter(exposures, Exposure > exposure_limit)

# detect signatures with no sample having at least the exposure_limit exposure to it, to show it in the plot
signatures_lower_exposure_than_limit = exposures %>% 
  group_by(Signature) %>% 
  summarise(max_exposure = max(Exposure)) %>% 
  filter(max_exposure < exposure_limit) %>% pull(Signature) %>% as.character

if(length(signatures_lower_exposure_than_limit) > 0){
  filtered_exposures = filtered_exposures %>% 
    bind_rows(data.frame(signatures_lower_exposure_than_limit) %>% 
                `colnames<-`("Signature"))}

## exposures (only samples with Exposures > exposure_limit%)
pos = position_jitter(w = 0.25, h = 0, seed = 1)
exposures_plot = ggplot(filtered_exposures, 
                        aes(x = Signature,
                            y = Exposure*100,
                            shape = dataset)) +
  scale_y_continuous(labels = function(x) sub("0+$", "", x)) +
  geom_point(aes(fill = info2),
             size = 4,
             position = pos) +
  scale_fill_manual(values = jet.colors(length(unique(exposures$info2)))) +
  scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(size=6, shape=21)),
         shape = guide_legend(override.aes = list(size=6))) +
  ggrepel::geom_text_repel(aes(label = paste(Sample, info1)),
                           size = 4,
                           force = 10,
                           position = pos,
                           max.overlaps = 1000000,
                           min.segment.length = 1) +
  facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
  theme_classic() +
  xlab("") +
  ggtitle(paste0("Model Score = ", round(good_model_score*100, 2), "%")) +
  ylab(paste0("% Exposure (>", exposure_limit*100, "%)")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(4, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

## weights
weights_plot = ggplot(weights %>%
                        mutate(`DNA repair\nactivity` = gsub("_2strands", "", `DNA repair\nactivity`),
                               Signature = factor(gsub("autoencoder", "", Signature), levels = seq(1:length(levels(weights$Signature))))) %>% 
                        # scale weights per signature so they add up to 1
                        group_by(Signature) %>% 
                        mutate(sumWeight = sum(Weight)) %>% 
                        group_by(Signature, `DNA repair\nactivity`) %>% 
                        summarise(Weight = Weight/sumWeight) %>% 
                        ungroup, 
                      aes(x = Signature,
                          y = Weight)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) sub("0+$", "", x)) +
  geom_col(aes(fill = `DNA repair\nactivity`)) +
  scale_fill_manual(values = jet.colors(length(levels(weights$`DNA repair\nactivity`)))) +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_grid(cols = vars(Signature), scales = "free") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2, "mm"),
        legend.text = element_text(size = 9))

combined_plots = cowplot::plot_grid(NULL,
                                    cowplot::plot_grid(exposures_plot, NULL, nrow = 1, rel_widths = c(1, 0.04*(optimal_k/11))),
                                    NULL,
                                    weights_plot,
                                    nrow = 4,
                                    rel_heights = c(0.02, 0.75,-0.05,1))
ggsave(paste0("autoencoder_exposures_weights_plot.jpg"),
       plot = combined_plots,
       device = "jpg",
       width = 21.3,
       height = 12,
       dpi = 700,
       bg = "white")
