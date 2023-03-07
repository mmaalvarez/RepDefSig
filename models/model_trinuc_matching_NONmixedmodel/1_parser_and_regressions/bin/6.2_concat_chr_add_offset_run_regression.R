### regression

# formula specifies which of the fold values that were used to simulate the mutcounts to use as dependent variable
formula = paste0(paste0("simulated_mutcount_", mutfoldinc, "fold ~ "),
                 paste(dnarep_marks$name, collapse = " + "), " + ",
                 "mb_domain + ",
                 "offset(log_freq_trinuc32)")

# stick to a generalized linear model for the negative binomial family
y = suppressWarnings(glm.nb(formula = formula, 
                            data = sim_pos_con))

## parse output
y_tidy = broom::tidy(y) %>%
  filter(term != "(Intercept)" & !str_detect(term, "mb_domain")) %>%
  # calc CI 95% from std error
  mutate("conf.low" = estimate - `std.error`*1.96,
         "conf.high" = estimate + `std.error`*1.96) %>% 
  select(c(term, estimate, contains("conf"))) %>% 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(sample_id = paste0(dnarep_mark_simulate, "-", altered_level, "__muts_", mutfoldinc, "-fold"),
         info1 = sample_id,
         info2 = sample_id,
         # theta value used, either because it was optimal or because it reached 10 iterations
         theta = y$theta) %>% 
  relocate(sample_id)
gc()

write_tsv(y_tidy, "simulated_positive_control.tsv")
