library(tidyverse)
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



jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))


##### parse input

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

# load results of regressions and just keep sample_id and coefficients for DNA repair marks
results_regressions = read_tsv("../../1_parser_and_regressions/res/results.tsv") %>%
  # make sure all samples are in the metadata table
  filter(sample_id %in% metadata$Sample) %>% 
  select(sample_id, contains("estimate_"), theta)



### PCA to look for biases regarding which samples were regressed with which theta values (either converged, or not)
pca_res = results_regressions %>% 
  select(sample_id, contains("estimate_")) %>% 
  column_to_rownames("sample_id") %>% 
  rename_with(~str_replace(., 'estimate_', '')) %>% 
  rename_with(~str_replace(., '.L', '')) %>%
  PCA(graph = FALSE)

pca = results_regressions %>%
  merge(pca_res$ind$coord %>% 
          data.frame %>% 
          rownames_to_column("sample_id")) %>% 
  select(sample_id, contains("Dim"), theta) %>% 
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  rename("Sample" = "sample_id") %>% 
  merge(metadata) %>% 
  ggplot(aes(x = PC1,
             y = PC2)) +
  coord_fixed() +
  geom_point(aes(fill = info2,
                 shape = as.character(floor(theta)))) +
  stat_ellipse(geom = "polygon",
               aes(col = as.character(floor(theta))),
               alpha = 0,
               show.legend = T,
               level = 0.95) +
  scale_fill_manual(values = jet.colors(length(unique(metadata$info2)))) +
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = jet.colors(length(unique(floor(results_regressions$theta))))) +
  guides(fill = guide_legend(override.aes = list(size=6, shape=21))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  ggrepel::geom_text_repel(aes(label = info1),
                           size = 1,
                           force = 4,
                           segment.size = 0.01,
                           min.segment.length = 0.001,
                           max.overlaps = 100000,
                           max.iter = 100000) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_blank())
ggsave("pca.svg",
       plot = pca,
       device = "svg",
       width = 20,
       height = 15.6,
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
ggsave("scree.svg",
       plot = scree,
       device = "svg",
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
ggsave("vars.svg",
       plot = vars,
       device = "svg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")
