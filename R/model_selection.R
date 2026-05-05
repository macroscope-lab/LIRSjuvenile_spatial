library(brms)
library(tidyverse)
library(tidybayes)
library(posterior)
library(loo)

# path to your saved models
model_dir <- "output/models"

# list all models
model_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)

# helper function to extract everything from one model
extract_model_info <- function(file) {
  
  fit <- readRDS(file)
  model_name <- basename(file) %>% str_remove(".rds")
  
  loo_obj <- loo(fit)
  
  loo_df <- tibble(
    model = model_name,
    elpd_loo = loo_obj$estimates["elpd_loo","Estimate"],
    p_loo    = loo_obj$estimates["p_loo","Estimate"],
    looic    = loo_obj$estimates["looic","Estimate"]
  )
  
  r2 <- bayes_R2(fit, summary = FALSE)
  
  r2_df <- tibble(
    model = model_name,
    R2_marginal = mean(r2),
    R2_marginal_low = quantile(r2, 0.025),
    R2_marginal_high = quantile(r2, 0.975)
  )
  
  # approximate conditional using re_formula = NA
  r2_cond <- bayes_R2(fit, re_formula = NA, summary = FALSE)
  
  r2_df <- r2_df %>%
    mutate(
      R2_conditional = mean(r2_cond),
      R2_cond_low = quantile(r2_cond, 0.025),
      R2_cond_high = quantile(r2_cond, 0.975)
    )

  fixef_draws <- as_draws_df(fit) %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value")
  
  fixef_summary <- fixef_draws %>%
    group_by(term) %>%
    summarise(
      mean = mean(value),
      
      l95 = quantile(value, 0.025),
      u95 = quantile(value, 0.975),
      
      l90 = quantile(value, 0.05),
      u90 = quantile(value, 0.95),
      
      l80 = quantile(value, 0.10),
      u80 = quantile(value, 0.90),
      
      .groups = "drop"
    ) %>%
    mutate(type = "fixed", model = model_name)

  ranef_draws <- as_draws_df(fit) %>%
    select(starts_with("r_")) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    # split term into components: group, level, coefficient
    extract(term, into = c("group", "level", "coef"),
            regex = "r_(.*?)\\[(.*?),(.*?)\\]")
  
  ranef_summary <- ranef_draws %>%
    group_by(group, level, coef) %>%
    summarise(
      mean = mean(value),
      
      l95 = quantile(value, 0.025),
      u95 = quantile(value, 0.975),
      
      l90 = quantile(value, 0.05),
      u90 = quantile(value, 0.95),
      
      l80 = quantile(value, 0.10),
      u80 = quantile(value, 0.90),
      
      .groups = "drop"
    ) %>%
    mutate(
      term = paste0("r_", group, "[", level, ",", coef, "]"),
      type = "random",
      model = model_name
    )
  
  effects_df <- bind_rows(fixef_summary, ranef_summary)
  
  # -------------------------
  # combine all summaries
  # -------------------------
  list(
    loo = loo_df,
    r2  = r2_df,
    effects = effects_df
  )
}

# -------------------------
# APPLY TO ALL MODELS
# -------------------------
results <- map(model_files, extract_model_info)

loo_table <- map_dfr(results, "loo")
r2_table  <- map_dfr(results, "r2")
effects_table <- map_dfr(results, "effects")

# -------------------------
# FINAL MODEL COMPARISON TABLE
# -------------------------
model_summary <- loo_table %>%
  left_join(r2_table, by = "model") %>%
  arrange(looic)

# view outputs
model_summary
effects_table

# optional: save
write_csv(model_summary, "output/dt1_cVSs_split_model_comparison_summary.csv")
write_csv(effects_table, "output/dt1_cVSs_split_model_effects_summary.csv")



library(tidyverse)

effects_clean <- effects_table %>%
  mutate(
    # clean names
    term = str_remove(term, "^b_"),
    term = str_replace_all(term, "_", " "),
    
    # identify random-effect groups
    group = case_when(
      type == "random" & str_detect(term, "site") ~ "site",
      type == "random" & str_detect(term, "year") ~ "year",
      TRUE ~ "fixed"
    )
  )

term_order <- effects_clean %>%
  group_by(term) %>%
  summarise(avg = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(avg) %>%
  pull(term)

effects_clean <- effects_clean %>%
  mutate(term = factor(term, levels = term_order))

x_limits <- range(
  effects_clean$l95,
  effects_clean$u95,
  na.rm = TRUE
)

library(ggplot2)

plot_effects <- function(df, effect_type) {
  
  df_sub <- df %>% filter(type == effect_type)
  
  p <- ggplot(df_sub, aes(x = mean, y = term)) +
    
    # 95% CI (thin)
    geom_errorbarh(aes(xmin = l95, xmax = u95),
                   height = 0, size = 0.4, alpha = 0.7) +
    
    # 80% CI (thick)
    geom_errorbarh(aes(xmin = l80, xmax = u80),
                   height = 0, size = 1.2) +
    
    geom_point(size = 1.8) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "grey40") +
    
    scale_x_continuous(limits = x_limits) +
    
    labs(x = NULL, y = NULL) +
    
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7)
    )
  
  # 👉 split random effects by group
  if (effect_type == "random") {
    p <- p + facet_wrap(~group, scales = "free_y")
  }
  
  return(p)
}

library(patchwork)

models <- unique(effects_clean$model)

plot_rows <- lapply(models, function(m) {
  
  df_m <- effects_clean %>% filter(model == m)
  
  p_fixed <- plot_effects(df_m, "fixed") +
    ggtitle(m) +
    theme(plot.title = element_text(size = 10, face = "bold"))
  
  p_random <- plot_effects(df_m, "random")
  
  p_fixed + p_random + plot_layout(widths = c(1, 1.2))
})

final_plot <- wrap_plots(plot_rows, ncol = 1)

final_plot


ggsave(
  "figures/dt1_cVSs_split_effects_caterpillar.png",
  final_plot,
  width = 12,
  height = 30,   # tall because lots of rows
  dpi = 300
)

