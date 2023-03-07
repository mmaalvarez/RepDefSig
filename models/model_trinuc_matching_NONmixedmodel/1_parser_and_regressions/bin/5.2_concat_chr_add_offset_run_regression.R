### regression

formula = paste0("mutcount ~ ",
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "mb_domain + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear model for the negative binomial family
y = suppressWarnings(glm.nb(formula = formula, 
                            data = sommut_dnarep_chromatin))

## parse output
y_tidy = broom::tidy(y) %>%
  filter(term != "(Intercept)" & !str_detect(term, "mb_domain")) %>%
  # calc CI 95% from std error
  mutate("conf.low" = estimate - `std.error`*1.96,
         "conf.high" = estimate + `std.error`*1.96) %>% 
  select(c(term, estimate, contains("conf"))) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(sample_id = sample,
         # theta value used, either because it was optimal or because it reached 10 iterations
         theta = y$theta)
gc()

## append features' coefficients and pvalues to metadata_sample
results_sample = full_join(metadata_sample, y_tidy) %>%
  relocate(sample_id) %>% 
  relocate(info1, info2, .before = theta)

write_tsv(results_sample, "results_sample.tsv")
