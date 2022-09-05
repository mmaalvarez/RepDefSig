library(tidyverse)
library(NMF) # for estimating optimal rank value (i.e. nº signatures)
library(RcppML) # for actual NMF -- dev version 0.5.4 (https://github.com/zdebruine/RcppML)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("nmf", "RcppML")


## parse input

metadata = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv") %>% 
  select(sample_id, source, tissue, OriginalType, hr_status, MSI_status, smoking_history, treatment_platinum, treatment_5FU, tumorPurity, gender) %>% 
  rename("Sample" = "sample_id")

# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../../1_parser_and_regressions/model_1/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_")) %>%  #  source, MSI_status, hr_status, smoking_history, treatment_platinum, treatment_5FU, tissue, OriginalType, 
  arrange(sample_id)

# a) keep positive coefficients, and convert negative to zero 
results_regressions_posmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))

# b) convert positive coefficients to zero, and convert negative to positive
results_regressions_negmatrix = results_regressions %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))

# merge them into the NMF input
coefficient_matrix = merge(results_regressions_posmatrix,
                           results_regressions_negmatrix,
                           by = "sample_id",
                           suffixes = c("_poscoeff", "_negcoeff")) %>%
  # transpose
  t() %>% 
  `colnames<-`(.[1, ]) %>% 
  .[-1, ] %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column("dna_repair_mark") %>% 
  mutate(dna_repair_mark = gsub("estimate_", "", dna_repair_mark)) %>% 
  column_to_rownames("dna_repair_mark") %>% 
  mutate_all(as.numeric) %>%
  as.matrix



### NMF

n_signatures = 10
# # get optimal n signatures ("rank")
# n_signatures = nmfEstimateRank(coefficient_matrix,
#                                range = 2,
#                                method = "brunet", # Brunet 2004 is the default --> nmf.getOption("default.algorithm")
#                                nrun = 2)
# NMF
nmf_res = nmf(coefficient_matrix,
              k = n_signatures, 
              seed = 1)

# Signature exposures in samples
exposures = nmf_res$h %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("Signature") %>%
  pivot_longer(cols = !contains("Signature"), names_to = "Sample", values_to = "Exposure") %>% 
  # add info on MSI, HR, smoking, and treatments
  left_join(metadata) %>% 
  mutate(Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h)))),
         MSI_parsed = ifelse(MSI_status %in% c("HYPER", "MSI", "ERCC2mut"), "MSI", NA),
         hr_parsed = ifelse(hr_status %in% c("HR_deficient"), "HRdef", NA),
         smoking_history_parsed = ifelse(smoking_history %in% c("Current", "Former"), "Smoker", NA),
         treatment_platinum_parsed = ifelse(treatment_platinum == "TRUE", "Platinum", NA),
         treatment_5FU_parsed = ifelse(treatment_5FU == "TRUE", "5FU", NA)) %>%
  unite(col = "Metadata", MSI_parsed, hr_parsed, smoking_history_parsed,treatment_platinum_parsed, treatment_5FU_parsed, na.rm = T, sep = " & ") %>% 
  mutate(Metadata = gsub("^$", NA, Metadata)) %>% 
  rename("Database" = "source") %>% 
  # highlight top hits for each signature
  group_by(Signature) %>% 
  mutate(is.hit = ifelse(Exposure==max(Exposure), "hit", NA),
         has.Metadata = ifelse(!is.na(Metadata), "yes", "no"))
write_tsv(exposures,
          "NMF_exposures.tsv")

# DNA repair mark weights in signatures
weights = nmf_res$w %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column("dna_repair_mark") %>%
  pivot_longer(cols = contains("nmf"), names_to = "signature", values_to = "weight") %>% 
  extract(dna_repair_mark, into = c("dna_repair_mark", "submatrix"), "(.*)_([^_]+$)") %>% 
  mutate(dna_repair_mark = factor(dna_repair_mark, levels = unique(gsub("_...coeff", "", rownames(nmf_res$w)))),
         signature = factor(signature, levels = unique(exposures$Signature))) %>% 
  arrange(dna_repair_mark, signature) %>% 
  group_by(dna_repair_mark, signature) %>% 
  ## keep larger of the 2 values for the 2 nmf submatrices (1 based on pos coefficients submatrix, and another on the absolute neg coefficients submatrix)
  summarise(Weight = max(weight)) %>% 
  ungroup %>% 
  rename("DNA repair activity" = "dna_repair_mark", 
         "Signature" = "signature") %>% 
  relocate(Signature)
write_tsv(weights,
          "NMF_weights.tsv")


### plotting

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## exposures (only Exposures > 0.001)
exposures_plot = ggplot(filter(exposures, Exposure > 0.001), 
                        aes(x = Signature,
                            # convert exposures to %
                            y = Exposure*100,
                            group = Database)) +
  scale_y_log10(labels = function(x) sub("0+$", "", x)) +
  # all points, no colors
  geom_point(aes(shape = Database),
             position = position_dodge(width = 1),
             alpha = 0.5,
             size = 4) +
  # colored those with metadata info
  geom_point(aes(shape = Database,
                 fill = Metadata,
                 alpha = has.Metadata),
             position = position_dodge(width = 1),
             size = 4) +
  scale_shape_manual(values = c(21,23,24,22,20,25)) + #"\u2716"
  scale_fill_manual(values = jet.colors(length(unique(filter(exposures, Exposure > 0.001 & !is.na(Metadata))$Metadata))),
                    na.value = "white") + # NA is assigned white
  scale_alpha_manual(values = c(0, 0.8), guide = 'none') +
  guides(shape = guide_legend(override.aes = list(size=6)),
         fill = guide_legend(override.aes = list(size=6, shape=21))) +
  # label sample names to top hits
  ggrepel::geom_text_repel(data = filter(exposures, is.hit == "hit"),
                           aes(label = Sample),
                           size = 3,
                           min.segment.length = 10000) +
  facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
  theme_classic() +
  xlab("") +
  ylab("% Exposure (>0.1% ; log10 scale)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1, "mm"))
  
## weights
weights_plot = ggplot(weights %>%
                        mutate(`DNA repair activity` = gsub("_2strands", "", `DNA repair activity`),
                               Signature = factor(gsub("nmf", "", Signature), levels = seq(1:length(levels(weights$Signature))))), 
                      aes(x = Signature,
                          y = Weight)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) sub("0+$", "", x)) +
  geom_col(aes(fill = `DNA repair activity`)) +
  scale_fill_manual(values = jet.colors(length(levels(weights$`DNA repair activity`)))) +
  guides(fill = guide_legend(override.aes = list(size=6))) +
  facet_grid(cols = vars(Signature), scales = "free", space = "free") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1, "mm"))

combined_plots = cowplot::plot_grid(NULL,
                                    exposures_plot,
                                    NULL,
                                    weights_plot,
                                    nrow = 4,
                                    rel_heights = c(0.02, 1,-0.05,1))
ggsave("NMF_plot.pdf",
       plot = combined_plots,
       device = "pdf",
       width = 22.5,
       height = 12,
       dpi = 600,
       bg = "white")

