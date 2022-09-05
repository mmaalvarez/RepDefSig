library(tidyverse)
library(NMF) # for estimating optimal rank value (i.e. nº signatures)
library(RcppML) # for actual NMF
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

# get optimal n signatures ("rank")
n_signatures = nmfEstimateRank(coefficient_matrix,
                               range = seq(1:100),
                               method = nmf.getOption("default.algorithm"),
                               nrun = 30)
# NMF
nmf_res = nmf(coefficient_matrix,
              k = n_signatures, 
              seed = 1)

# Signature exposures in samples
exposures = nmf_res$h %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("Signature") %>%
  pivot_longer(cols = !contains("Signature"), names_to = "Sample", values_to = "Exposure") %>% 
  left_join(metadata) %>% 
  mutate(Signature = factor(Signature, levels = unique(exposures$Signature)))

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


### plotting

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## exposures (only Exposures > 0.001)
exposures_plot = ggplot(filter(exposures, Exposure>0.001), 
                        aes(x = Signature,
                            y = Exposure,
                            shape = source)) +
  scale_y_sqrt(breaks = seq(0, 1, 0.001)) +
  ggbeeswarm::geom_beeswarm(alpha = 0.5,
                            dodge.width = 1) +
  ggrepel::geom_text_repel(data = filter(exposures, Exposure>0.004),
                           aes(label = Sample),
                           force = 10,
                           size = 1,
                           max.overlaps = 1000) +
  facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
  theme_classic() +
  xlab("") +
  ylab("Sample exposures >0.001 (sqrt. scale)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_blank())
  
## weights
weights_plot = ggplot(weights, 
                      aes(x = Signature,
                          y = Weight,
                          fill = `DNA repair activity`)) +
  geom_col() +
  scale_fill_manual(values = jet.colors(length(levels(weights$`DNA repair activity`)))) +
  facet_grid(cols = vars(Signature), scales = "free", space = "free") +
  theme_classic() +
  ylab("DNA repair activity weights") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_blank())

combined_plots = cowplot::plot_grid(exposures_plot,
                                    weights_plot,
                                    nrow = 2)
ggsave("model1_NMF.pdf",
       plot = combined_plots,
       device = "pdf",
       width = 22,
       height = 12,
       dpi = 600,
       bg = "white")
