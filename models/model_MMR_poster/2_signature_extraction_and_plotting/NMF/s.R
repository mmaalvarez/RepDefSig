library(tidyverse)
library(RcppML) # for NMF
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


dir.create("plots")


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


## CTCF + cohesin
CTCFcohesin_CC2TT_genomewide = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/2_count_CC2TT_genomewide_in_tumors/res/samples_CC2TT_genomewide_counted.tsv")
  
## N A3A-putative TpC2DpH mut pairs IN SKIN (to compare to CTCF's CC>TT caused by UV in skin, but A3A can be confounded by UV, i.e. in skin)
TpC2DpH_mutpairs_SKIN_SKIN = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/A3A_hairpins/2_SKIN_count_1kbpairs_TpC2DpH_in_tumors/res/samples_A3A_1kbpairs_TpC2DpH_counted.tsv")


metadata = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/RPE1_POLH_WTvsKO/1_parse_vcfs/sample_treatments.tsv",
             "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/metadatacomb_metadata_final_6datasets__noconsent_44plus11_samples_removed.csv") %>% 
  # only sample_id and info* columns are selected
  map_df(~read_tsv(.x)) %>% 
  mutate(dataset = ifelse(str_detect(sample_id, "MSM0"),
                          "Kucab et al. 2019",
                          ifelse(str_detect(sample_id, "MSK0"),
                                 "Zou et al. 2021",
                                 ifelse(str_detect(sample_id, "POLH"),
                                        "RPE1 POLHwt/ko ctrl/UV",
                                        source))),
         info1 = ifelse(str_detect(sample_id, "MSK0"),
                        paste0(info1, "ko"),
                        info1),
         info2 = gsub("^[a-z]_", "", info2), 
         info2 = ifelse(str_detect(sample_id, "MSM0") | str_detect(sample_id, "MSK0"),
                        gsub("DNA damage response inhibitors", "DNA damage resp. inh.", info2),
                        info2),
         info2 = ifelse(str_detect(sample_id, "MSK0"),
                        paste0(info2, " (pathway)"),
                        ifelse(str_detect(sample_id, "MSM0"),
                               paste0(info2, " (treatment)"),
                               info2))) %>% 
  rename("sample_id2" = "sample_id") %>% 
  merge(names_conversion_table, all = T) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            ifelse(!is.na(sample_id2),
                                   sample_id2,
                                   sample_id_2),
                            sample_id)) 

metadata_CTCFcohesin_CC2TT_genomewide = metadata %>% 
  merge(CTCFcohesin_CC2TT_genomewide, all = T) %>% 
  select(-c(sample_id, sample_id_2)) %>% 
  rename("sample_id" = "sample_id2") %>% 
  merge(CTCFcohesin_CC2TT_genomewide, all = T) %>% 
  merge(names_conversion_table, all = T) %>% 
  select(-sample_id2) %>% 
  distinct
  
metadata_TpC2DpH_mutpairs_SKIN = metadata %>% 
  merge(TpC2DpH_mutpairs_SKIN, all = T) %>% 
  mutate(sample_id = ifelse(sample_id!=sample_id2 & !is.na(sample_id2),
                            sample_id2,
                            sample_id)) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarise_all(funs(toString(na.omit(.)))) %>% 
  mutate_all(na_if, "") %>% 
  distinct

metadata_CC2TT_TpC2DpH = merge(metadata_CTCFcohesin_CC2TT_genomewide,
                               metadata_TpC2DpH_mutpairs_SKIN, all = T) %>% 
  filter(!(is.na(CC2TT_genomewide) & is.na(A3A_1kbpairs_TpC2DpH))) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarise_all(funs(toString(na.omit(.)))) %>% 
  mutate_all(na_if, "")



## results regressions

results_regressions = lapply(c("../../1_parser_and_regressions/res/results_real_samples.tsv",
                               "../../1_parser_and_regressions/res/simulated_positive_controls.tsv"),
                             read_tsv) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>% 
  rename_with(~str_replace(., 'high$', '')) %>% 
  rename_with(~str_replace(., 'hairpin_TpCpH$', '')) %>% 
  rename_with(~str_replace(., 'CTCF_cohesin_peak$', '')) %>%
  rename("sample_id2" = "sample_id") %>% 
  left_join(metadata_CC2TT_TpC2DpH) %>% 
  mutate(sample_id = ifelse(is.na(sample_id),
                            sample_id2,
                            sample_id)) %>% 
  select(-c(info1,info2)) %>% 
  drop_na(starts_with("estimate_")) %>% 
  select(sample_id, contains("estimate_"), contains("conf"))

coefficient_table = results_regressions %>%
  rename_all(~str_replace_all(., 'estimate_', 'estimate ')) %>% 
  rename_all(~str_replace_all(., 'conf.high_', 'conf.high ')) %>% 
  rename_all(~str_replace_all(., 'conf.low_', 'conf.low ')) %>% 
  pivot_longer(cols = -sample_id , names_to = 'stat_mark', values_to = 'value') %>% 
  separate(stat_mark, into = c("stat", "mark"), sep = " ") %>% 
  arrange(sample_id, mark) %>% 
  pivot_wider(names_from = stat) %>% 
  group_by(sample_id, mark)


## Parameters and initializing of some objects
set.seed(1)
maxK = length(unique(coefficient_table$mark)) # max number of signatures to consider, it will go from 2 to maxK -- shouldn't be larger than nº of features



#################################################################################
###### NMF


##### split the original coefficients between pos and neg
# a) keep positive coefficients, and convert negative to zero 
coefficients_posmatrix = coefficient_table %>% 
  select(sample_id, mark, contains("estimate")) %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))
# b) convert positive coefficients to zero, and convert negative to positive
coefficients_negmatrix = coefficient_table %>% 
  select(sample_id, mark, contains("estimate")) %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))

# regenerate the coeff matrix and transpose, for RcppML::nmf()
coefficient_matrix_RcppML = bind_rows(mutate(coefficients_posmatrix, submatrix = "poscoeff"),
                                      mutate(coefficients_negmatrix, submatrix = "negcoeff")) %>% 
  unite("dna_repair_mark",mark, submatrix, sep = "_") %>% 
  mutate(dna_repair_mark = gsub("estimate_", "", dna_repair_mark)) %>% 
  # transpose
  pivot_wider(names_from = sample_id, values_from = estimate) %>% 
  column_to_rownames("dna_repair_mark") %>% 
  as.matrix


##### prepare condition_pathway_pairs for "good-model-score"
repair_mark_pathways = coefficient_table %>% 
  select(sample_id,  mark, estimate) %>% 
  rename("dna_repair_mark" = "mark",
         "sample_id" = "sample_id") %>% 
  mutate(putative_repair_pathway_involved = ifelse(str_detect(dna_repair_mark, "OGG1_"),
                                                   "BER",
                                                   ifelse(str_detect(dna_repair_mark, "UV_"),
                                                          "NER",
                                                          ifelse(str_detect(dna_repair_mark, "MSH6_|SETD2_"),
                                                                 "MMR",
                                                                 ifelse(str_detect(dna_repair_mark, "XRCC4"),
                                                                        "DSBR",
                                                                        "other"))))) # ,#ADD TP53
condition_pathway_pairs = metadata_CC2TT_TpC2DpH %>% 
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
                                                                               "other")))))) %>% # ADD TP53
  separate_rows(putative_repair_pathway_involved, sep = ",") %>% 
  filter(!is.na(putative_repair_pathway_involved)) %>% 
  left_join(metadata_CC2TT_TpC2DpH) %>% 
  filter(sample_id %in% coefficient_table$sample_id) %>% 
  arrange(putative_repair_pathway_involved) %>% 
  left_join(repair_mark_pathways) %>% 
  select(sample_id, dna_repair_mark) %>% 
  rename("DNA repair\nactivity" = "dna_repair_mark") %>% 
  group_by(sample_id) %>%
  mutate(sample_id_possible_marks = n()) %>%
  ungroup()
condition_pathway_pairs = condition_pathway_pairs %>% 
  rowwise %>% 
  mutate(max_sample_mark_score = 1 / length(unique(condition_pathway_pairs$sample_id)) / sample_id_possible_marks) %>%
  select(-c(sample_id_possible_marks))


####
### Run NMF for several k's

# for plots
jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))

for(optimal_k in seq(3, maxK)){
  
  # final NMF (here using RcppML::nmf instead of NMF::nmf as above)
  nmf_res = RcppML::nmf(coefficient_matrix_RcppML, 
                        k = optimal_k,
                        maxit = 10000, 
                        seed = 1)
  
  # Signature exposures in samples
  exposures = nmf_res$h %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column("Signature") %>%
    pivot_longer(cols = !contains("Signature"), names_to = "sample_id", values_to = "Exposure") %>% 
    # add metadata info (e.g. treatments, MSI, HR, smoking...)
    merge(metadata_CC2TT_TpC2DpH, all = T) %>% 
    filter(!is.na(Exposure)) %>% 
    separate(sample_id, into = c("sample_id", "mut x-fold"), sep = "__", fill = "right") %>% 
    # simulated pos controls dont have CC>TT nor A3A mut pairs, so add them an average fake one so that their size is not tiny in the plots
    mutate(CC2TT_genomewide = ifelse(is.na(CC2TT_genomewide) & !is.na(`mut x-fold`),
                                     mean(metadata$CC2TT_genomewide, na.rm = T),
                                     CC2TT_genomewide)) %>%
    mutate(A3A_1kbpairs_TpC2DpH = ifelse(is.na(A3A_1kbpairs_TpC2DpH) & !is.na(`mut x-fold`),
                                        mean(metadata$A3A_1kbpairs_TpC2DpH, na.rm = T),
                                        A3A_1kbpairs_TpC2DpH)) %>%
    mutate(sample_id = gsub("-high_activity", "", sample_id),
           dataset = ifelse(is.na(dataset),
                            `mut x-fold`,
                            dataset),
           Signature = factor(Signature, levels = unique(rownames(data.frame(nmf_res$h))))) %>%
    # highlight top hits for each signature
    group_by(Signature) %>% 
    mutate(is.hit = ifelse(Exposure==max(Exposure), "Top hit", NA))
  
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
    ## subtract results of pos - results of neg, and get the absolute
    summarise(Weight = abs(.Primitive("-")(weight[1], weight[2]))) %>% 
    ungroup %>% 
    rename("DNA repair\nactivity" = "dna_repair_mark", 
           "Signature" = "signature") %>% 
    relocate(Signature)
  
  #####
  ## calculate "good-model-score": 1 would mean that ALL samples have 100% exposure from the signature(s) that is contributed 100% by a specific DNArep mark whose mechanism/pathway overlaps completely with the sample´s condition (gene-/-, treatment...), and 0 the opposite
  # actually only considering "ground-truth" sample_condition--DNArep_pathway pairs ('condition_pathway_pairs')
  good_model_score_table = weights %>% 
    merge(condition_pathway_pairs, all = T) %>% 
    merge(exposures, all = T) %>% 
    rename("sensical sample-mark pair" = "DNA repair\nactivity",
           "mark weight" = "Weight",
           "signature exposure" = "Exposure") %>% 
    select(sample_id, "sensical sample-mark pair", Signature, "mark weight", "signature exposure", max_sample_mark_score) %>% 
    mutate(sample_mark_score = `mark weight` * `signature exposure` * max_sample_mark_score) %>% 
    arrange(sample_id, `sensical sample-mark pair`, Signature)
  
  #write_tsv(good_model_score_table, paste0("K", optimal_k, "_table.tsv"))
  
  good_model_score = drop_na(good_model_score_table) %>% pull(sample_mark_score) %>% sum
  #####
  
  
  ### plotting
  
  # for each nmf signature, plot only the top n samples with highest Exposure
  top_n_samples = 5
  filtered_exposures = exposures %>% 
    group_by(Signature) %>% 
    arrange(desc(Exposure)) %>% 
    slice_head(n = top_n_samples) %>% 
    ungroup
  
  ## exposures (only top_n_samples)
  pos = position_jitter(w = 0.25, h = 0, seed = 1)
  exposures_plot = ggplot(filtered_exposures %>% 
                            mutate(`Genome-wide CC>TT` = as.numeric(CC2TT_genomewide),
                                   `5'-T(C>D)H-3' pairs` = as.numeric(A3A_1kbpairs_TpC2DpH)), 
                          aes(x = Signature,
                              y = Exposure*100)) +
    scale_y_continuous(labels = function(x) sub("0+$", "", x)) +
    geom_point(aes(fill = `Genome-wide CC>TT`,
                   size = `5'-T(C>D)H-3' pairs`),
                   shape = 21, #dataset, ##scale_shape_manual(values = c(21, 23, 25, seq(1,length(unique(filtered_exposures$dataset))-3,1))) +
                   position = pos) +
    scale_fill_gradient(low = "white", high = "black") +
    guides(fill = guide_legend(override.aes = list(size=4, shape=21)),
           shape = guide_legend(override.aes = list(size=4))) +
    ggrepel::geom_text_repel(aes(label = sample_id),
                             size = 3,
                             force = 5,
                             position = pos,
                             max.overlaps = 1000000,
                             min.segment.length = 1) +
    facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
    theme_classic() +
    xlab("") +
    #ggtitle(paste0("Model Score = ", round(good_model_score*100, 2), "%")) +
    ylab(paste0("% Exposure (top-", top_n_samples, " samples)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(6, "mm"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 8))
  
  ## weights
  weights_plot = ggplot(weights %>%
                          mutate(`DNA repair\nactivity` = gsub("_2strands", "", `DNA repair\nactivity`),
                                 Signature = factor(gsub("nmf", "", Signature), levels = seq(1:length(levels(weights$Signature))))) %>% 
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
          legend.text = element_text(size = 8))
  
  combined_plots = cowplot::plot_grid(NULL,
                                      cowplot::plot_grid(exposures_plot, NULL, nrow = 1, rel_widths = c(1, 0.04*(optimal_k/11))),
                                      NULL,
                                      weights_plot,
                                      nrow = 4,
                                      rel_heights = c(0.02, 0.75,-0.05,1))
  ggsave(paste0("plots/NMF_exposures_weights_plot_k", optimal_k, ".jpg"),
         plot = combined_plots,
         device = "jpg",
         width = 21.3,
         height = 12,
         dpi = 700,
         bg = "white")
}
