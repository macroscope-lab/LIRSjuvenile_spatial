library(brms)
library(tidyverse)
library(tidybayes)
library(posterior)
library(loo)
library(ggplot2)
library(patchwork)
library(rstan)
## helper function to extract everything from one model ####
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
  
  # Bayesian R2 analogues for brms mixed models
  # Marginal R2 analogue: fixed effects only
  r2_marg <- as.numeric(bayes_R2(
    fit,
    re_formula = NA,
    summary = FALSE
  ))
  
  # Conditional R2 analogue: fixed + random effects
  r2_cond <- as.numeric(bayes_R2(
    fit,
    re_formula = NULL,
    summary = FALSE
  ))
  
  r2_df <- tibble(
    model = model_name,
    
    R2_marginal = mean(r2_marg),
    R2_marginal_low = quantile(r2_marg, 0.025),
    R2_marginal_high = quantile(r2_marg, 0.975),
    
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
    tidyr::extract(term, into = c("group", "level", "coef"),
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
  
  # combine all summaries
  list(
    loo = loo_df,
    r2  = r2_df,
    effects = effects_df
  )
}

# path to your saved models ####
setwd("~/Documents/LIRSjuvenile_spatial")
model_dir <- "output/models"
# list all models
model_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
model_files <- model_files[grepl("aint_ry", model_files)]
model_files <- model_files[-3]

results <- map(model_files, extract_model_info)

loo_table <- map_dfr(results, "loo")
r2_table  <- map_dfr(results, "r2")
effects_table <- map_dfr(results, "effects")

# FINAL MODEL COMPARISON TABLE
model_summary <- loo_table %>%
  left_join(r2_table, by = "model") %>%
  arrange(looic)

# view outputs
model_summary
effects_table

# save
# write_csv(model_summary, "output/intVSrandomVBadu_model_comparison_summary.csv")
#write_csv(effects_table, "output/dt1_cVSs_split_model_effects_summary.csv")

effects_clean <- effects_table %>%
  filter(term != "b_Intercept")%>%
  mutate(
    scale = str_sub(model,start=10 , end=10),
    # clean names
    term = str_remove(term, "^b_"),
    term = str_replace_all(term, "_", " "),
    scale = str_replace_all(scale, "1", "scale = 1m x 1m"),
    scale = str_replace_all(scale, "2", "scale = 2m x 2m"),
    scale = str_replace_all(scale, "5", "scale = 0.5m x 0.5m"),
    # identify random-effect groups
    group = case_when(
      type == "random" & str_detect(term, "site") ~ "site",
      type == "random" & str_detect(term, "year") ~ "year",
      TRUE ~ "fixed"
    ),
    var = case_when(
      str_detect(model, "full") ~ "full model",
      TRUE ~ "alone"
    ),
    # more cleaning
    term = str_remove(term, "r site\\["),
    term = str_remove(term, "r year\\["),
    term = str_remove(term, ",Intercept\\]"),   
    term = str_replace_all(term, "Dc", "fractal dimension"),
    term = str_replace_all(term, "R log10 c","surface rugosity (log₁₀)"),
    term = str_replace_all(term, "mean sett log c", "mean settlement (log)"),
    term = str_replace_all(term, "n adu log c", "adult abundance (log)"),
    
    
  )

term_order <- effects_clean %>%
  group_by(term) %>%
  summarise(avg = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(avg) %>%
  pull(term)

effects_clean <- effects_clean %>%
  mutate(term = factor(term, levels = term_order))

x_limits <- range(
  effects_clean[effects_clean$type == "fixed",]$l95,
  effects_clean[effects_clean$type == "fixed",]$u95,
  na.rm = TRUE
)

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
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 6)
    )
  
  # split random effects by group
  if (effect_type == "random") {
    p <- p + facet_wrap(~group, scales = "free_y")
  }
  
  return(p)
}

models <- unique(effects_clean$model)

title_df <- tibble::tribble(
  ~model, ~title,
  "brm_nb_dt5_full_c_aint_ry", "0.5x0.5 m scale",
  "brm_nb_dt1_full_c_aint_ry", "1x1m scale",
  "brm_nb_dt2_full_c_aint_ry", "2x2m scale",
  "brm_nb_dt4_full_c_aint_ry", "4x4m scale"
)

plot_rows <- lapply(models, function(m) {
  
  df_m <- effects_clean %>%
    filter(model == m)
  
  plot_title <- title_df %>%
    filter(model == m) %>%
    pull(title)
  
  p_fixed <- plot_effects(df_m, "fixed") +
    ggtitle(plot_title) +
    scale_x_continuous(limits = NULL) +
    theme(plot.title = element_text(
      size = 8,
      face = "bold")
    )
  
  p_random <- plot_effects(df_m, "random")
  
  p_fixed + p_random +
    plot_layout(widths = c(1, 1.2))
})

final_plot <- wrap_plots(plot_rows, ncol = 1)
final_plot

ggsave("figures/models_effects_caterpillar_across_scales.pdf",
       final_plot,width = 5,height = 5,dpi = 300)
ggsave("figures/models_effects_caterpillar_across_scales.png",
       final_plot,width = 5,height = 5,dpi = 300)

## effect figure #####

plot_effect <- function(mod_file,
                        dt_file,
                        predictor_raw,
                        predictor_c,
                        xlab,
                        response = "n_juv") {
  
  mod <- mod_file
  dt  <- dt_file
  
  ## Store raw means used for centering
  n_adu_log_mean     <- mean(dt$n_adu_log, na.rm = TRUE)
  mean_sett_log_mean <- mean(dt$mean_sett_log, na.rm = TRUE)
  R_log10_mean       <- mean(dt$R_log10, na.rm = TRUE)
  D_mean             <- mean(dt$D, na.rm = TRUE)
  
  ## Create centered variables
  dt$n_adu_log_c     <- dt$n_adu_log - n_adu_log_mean
  dt$mean_sett_log_c <- dt$mean_sett_log - mean_sett_log_mean
  dt$R_log10_c       <- dt$R_log10 - R_log10_mean
  dt$Dc              <- dt$D - D_mean
  
  ## Build prediction grid on the original/raw scale
  x_raw <- seq(
    min(dt[[predictor_raw]], na.rm = TRUE),
    max(dt[[predictor_raw]], na.rm = TRUE),
    length.out = 100
  )
  
  ## Prediction data: model still receives centered predictors
  newdata <- tibble(
    n_adu_log_c      = 0,
    mean_sett_log_c  = 0,
    R_log10_c        = 0,
    Dc               = 0,
    x_raw_plot       = x_raw
  )
  
  ## Convert the focal raw predictor into the centered version
  if (predictor_c == "n_adu_log_c") {
    newdata[[predictor_c]] <- x_raw - n_adu_log_mean
  }
  
  if (predictor_c == "mean_sett_log_c") {
    newdata[[predictor_c]] <- x_raw - mean_sett_log_mean
  }
  
  if (predictor_c == "R_log10_c") {
    newdata[[predictor_c]] <- x_raw - R_log10_mean
  }
  
  if (predictor_c == "Dc") {
    newdata[[predictor_c]] <- x_raw - D_mean
  }
  
  ## Posterior predictions
  preds <- newdata |>
    add_epred_draws(
      mod,
      re_formula = NA,
      ndraws = 1000
    ) |>
    median_qi(
      .epred,
      .width = c(.5, .8, .95)
    )
  
  ## Plot
  ggplot(
    preds,
    aes(
      x = x_raw_plot,
      y = .epred,
      ymin = .lower,
      ymax = .upper
    )
  ) +
    
    geom_point(
      data = dt,
      aes(
        x = .data[[predictor_raw]],
        y = .data[[response]]
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      alpha = .20,
      size = 1.2
    ) +
    
    geom_lineribbon(
      aes(fill = factor(.width)),
      alpha = .50
    ) +
    
    geom_line(
      colour = "#3A2767",
      linewidth = .8
    ) +
    
    scale_fill_manual(
      values = c(
        "0.5" = "#3A2767",
        "0.8" = "#7A5BA8",
        "0.95" = "#B19CD9"
      ),
      name = "Credible interval"
    ) +
    
    labs(
      x = xlab,
      y = "Predicted juvenile abundance"
    ) +
    
    theme_classic(base_size = 14) +
    
    coord_cartesian(ylim = c(0, 110)) +
    
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
}


scales <- data.frame(
  mod = c( "output/models/brm_nb_dt1_full_c_aint_ry.rds","output/models/brm_nb_dt2_full_c_aint_ry.rds","output/models/brm_nb_dt5_full_c_aint_ry.rds"),
  dt = c("output/dt1.csv","output/dt2.csv","output/dt50.csv"),
  out = c("figures/fef_dt1_aint_ry.pdf","figures/fef_dt2_aint_ry.pdf","figures/fef_dt50_aint_ry.pdf"),
  out_png = c("figures/fef_dt1_aint_ry.png","figures/fef_dt2_aint_ry.png","figures/fef_dt50_aint_ry.png"))

for(i in seq_len(nrow(scales))) {
  mod <- readRDS(scales$mod[i])
  dt  <- read.csv(scales$dt[i])
  
  p1 <- plot_effect(
    mod,
    dt,
    predictor_raw = "n_adu_log",
    predictor_c = "n_adu_log_c",
    xlab = "adult abundance (log)"
  )
  
  p2 <- plot_effect(
    mod,
    dt,
    predictor_raw = "n_adu_log",
    predictor_c = "n_adu_log_c",
    xlab = "mean settlement (log)"
  )
  
  
  p3 <- plot_effect(
    mod,
    dt,
    predictor_raw = "n_adu_log",
    predictor_c = "n_adu_log_c",
    xlab = "surface rugosity (log₁₀)"
  )
  
  p4 <- plot_effect(
    mod,
    dt,
    predictor_raw = "n_adu_log",
    predictor_c = "n_adu_log_c",
    xlab = "fractal dimension"
  )
  
  fef <- (p1|p2)/(p3|p4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  fef
  ggsave(scales$out[i], fef, width = 8, height = 8)
  ggsave(scales$out_png[i], fef, width = 8, height = 8)
  
  
}


### interaction raster plots with uncertainty #####

make_interaction_plot <- function(mod_file, dt_file, outfile, outfile_png) {
  
  library(brms)
  library(ggplot2)
  library(patchwork)
  
  mod <- readRDS(mod_file)
  dt  <- read.csv(dt_file)
  dt$site <- factor(dt$site)
  dt$year <- factor(dt$year)
  dt$n_adu_log_c     <- as.numeric(scale(dt$n_adu_log,     scale = FALSE))
  dt$mean_sett_log_c <- as.numeric(scale(dt$mean_sett_log, scale = FALSE))
  dt$R_log10_c       <- as.numeric(scale(dt$R_log10,       scale = FALSE))
  dt$Dc              <- as.numeric(scale(dt$D,             scale = FALSE))
  
  preds <- c("n_adu_log_c", "Dc", "R_log10_c", "mean_sett_log_c")
  
  labs_map <- c(
    n_adu_log_c     = "adult abundance (log, centered)",
    Dc              = "fractal dimension (centered)",
    R_log10_c       = "surface rugosity (log₁₀, centered)",
    mean_sett_log_c = "average settlers (log, centered)"
  )
  
  # build aggregated prediction data for each pair
  build_pair_data <- function(var_x, var_y, n_grid = 80) {
    
    grid_x <- seq(min(dt[[var_x]], na.rm = TRUE),
                  max(dt[[var_x]], na.rm = TRUE),
                  length.out = n_grid)
    grid_y <- seq(min(dt[[var_y]], na.rm = TRUE),
                  max(dt[[var_y]], na.rm = TRUE),
                  length.out = n_grid)
    
    pred_grid <- expand.grid(
      year  = levels(dt$year),
      x_val = grid_x,
      y_val = grid_y
    )
    
    for (p in preds) pred_grid[[p]] <- 0
    pred_grid[[var_x]] <- pred_grid$x_val
    pred_grid[[var_y]] <- pred_grid$y_val
    
    fit <- fitted(mod, newdata = pred_grid, re_formula = NA, summary = TRUE)
    pred_grid$estimate <- fit[, "Estimate"]
    pred_grid$ci_width <- fit[, "Q97.5"] - fit[, "Q2.5"]
    
    agg <- aggregate(cbind(estimate, ci_width) ~ get(var_x) + get(var_y),
                     data = pred_grid, FUN = mean)
    names(agg)[1:2] <- c(var_x, var_y)
    agg
  }
  
  other_preds <- setdiff(preds, "n_adu_log_c")
  pair_data <- lapply(other_preds, function(v) {
    list(var_y = v, data = build_pair_data("n_adu_log_c", v))
  })
  
  # compute shared color scale limits?
  est_range <- range(unlist(lapply(pair_data, function(x) x$data$estimate)),
                     na.rm = TRUE)
  ci_range  <- range(unlist(lapply(pair_data, function(x) x$data$ci_width)),
                     na.rm = TRUE)
  
  # build panels using shared color scale
  make_panel <- function(agg, var_y, fill_var, fill_lab, palette, limits, title) {
    ggplot(agg, aes_string(x = "n_adu_log_c", y = var_y, fill = fill_var)) +
      geom_raster(interpolate = TRUE, alpha = .95) +
      geom_contour(aes_string(z = fill_var),
                   colour = "white", alpha = 0.4, linewidth = 0.3) +
      geom_point(data = dt,
                 aes_string(x = "n_adu_log_c", y = var_y),
                 inherit.aes = FALSE,
                 colour = "white", alpha = 0.2, size = 1) +
      scale_fill_viridis_c(option = palette, name = fill_lab, limits = limits) +
      labs(x = labs_map[["n_adu_log_c"]], y = labs_map[[var_y]],
           title = title) +
      theme_classic(base_size = 13) +
      theme(
        axis.text    = element_text(colour = "black"),
        panel.border = element_blank()
      )
  }
  
  pair_rows <- lapply(pair_data, function(pd) {
    mean_panel <- make_panel(pd$data, pd$var_y,
                             "estimate", "Juvenile\nabundance",
                             "magma", est_range,
                             "Posterior mean")
    ci_panel   <- make_panel(pd$data, pd$var_y,
                             "ci_width", "95% CI width",
                             "viridis", ci_range,
                             "Uncertainty")
    mean_panel + ci_panel
  })
  
  # stack rows and collect shared legends
  p <- wrap_plots(pair_rows, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Interactions with adult abundance",
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    ) &
    theme(legend.position = "right")
  
  ggsave(outfile, p, width = 11, height = 13)
  ggsave(outfile_png, p, width = 11, height = 13)
  return(p)
}

scales <- data.frame(
  mod = c( "output/models/brm_nb_dt1_full_c_aint_ry.rds","output/models/brm_nb_dt2_full_c_aint_ry.rds","output/models/brm_nb_dt5_full_c_aint_ry.rds"),
  dt = c("output/dt1.csv","output/dt2.csv","output/dt50.csv"),
  out = c("figures/interaction_panels_dt1_aint_ry.pdf","figures/interaction_panels_dt2_aint_ry.pdf","figures/interaction_panels_dt50_aint_ry.pdf"),
  out_png = c("figures/interaction_panels_dt1_aint_ry.png","figures/interaction_panels_dt2_aint_ry.png","figures/interaction_panels_dt50_aint_ry.png"),
  out_pd = c("figures/effSizes_and_postDraws_brm_nb_dt1_full_c_aint_ry.pdf", "figures/effSizes_and_postDraws_brm_nb_dt2_full_c_aint_ry.pdf",
            "figures/effSizes_and_postDraws_brm_nb_dt5_full_c_aint_ry.pdf"))

preds <- c("n_adu_log_c","mean_sett_log_c","R_log10_c","Dc")

for(i in seq_len(nrow(scales))) {
  make_interaction_plot(
    mod_file = scales$mod[i],
    dt_file  = scales$dt[i],
    outfile  = scales$out[i],
    outfile_png = scales$out_png[i]
  )
}

### make interaction and effects plots ####
  mod_file <-  scales$mod[i]
  dt_file <- scales$dt[i]
  out_file <- scales$out_pd[i]

  df_m <- effects_clean %>%
    filter(model == substr(mod_file,15,39)) %>%
    mutate(
      term_label = recode(
        term,
        "adults (log)" = "adult abundance",
        "settlers (log)" = "average settlers",
        "R log10" = "surface rugosity",
        "D" = "fractal dimension",
        "year2023" = "Year (2023)",
        "adults (log):settlers (log)" = "adult abundance × average settlers",
        "adults (log):R log10" = "adult abundance × surface rugosity",
        "adults (log):D" = "adult abundance × fractal dimension",
      ),
      effect_type = if_else(
        grepl(":", term),
        "Interaction",
        "Main effect"
      )
    )
  
  df_m <- df_m %>%
    arrange(desc(effect_type == "Interaction"),desc(mean)) %>%
    mutate(term_label = factor(term_label,levels = rev(term_label))) %>%
    filter(type == "fixed")
  e <- ggplot(df_m, aes(x = mean, y = term_label)) +
    geom_errorbarh(
      aes(xmin = l95, xmax = u95),
      height = 0,
      size = 0.5,
      alpha = 0.7
    ) +
    geom_errorbarh(
      aes(xmin = l80, xmax = u80),
      height = 0,
      size = 1.5
    ) +
    geom_point(size = 3.5) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      size = 1.2
    ) +
    scale_x_continuous() +
    labs(
      x = "Model effect sizes",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 11),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 14),
      plot.margin = margin(t = 12, r = 8, b = 8, l = 8)
    )
  
  
  mod <- readRDS(mod_file)
  draws <- as_draws_df(mod)
  
  draws_long <- tibble(
    # Main effects
    `adult abundance`                     = draws$b_n_adu_log_c,
    `average settlers`                    = draws$b_mean_sett_log_c,
    `surface rugosity`                    = draws$b_R_log10_c,
    `fractal dimension`                   = draws$b_Dc,
    
    # Interactions with adult abundance
    `adult abundance x average settlers`  = draws$`b_n_adu_log_c:mean_sett_log_c`,
    `adult abundance x surface rugosity`  = draws$`b_n_adu_log_c:R_log10_c`,
    `adult abundance x fractal dimension` = draws$`b_n_adu_log_c:Dc`
  ) %>%
    pivot_longer(
      everything(),
      names_to = "predictor",
      values_to = "estimate"
    ) %>%
    mutate(
      type = if_else(str_detect(predictor, " x "), "Interaction", "Main effect"),
      predictor = factor(
        predictor,
        levels = c(
          "adult abundance",
          "average settlers",
          "surface rugosity",
          "fractal dimension",
          "adult abundance x average settlers",
          "adult abundance x surface rugosity",
          "adult abundance x fractal dimension"
        )
      )
    )
  
  
  d <- ggplot(draws_long, aes(x = estimate, y = fct_rev(predictor), fill = type)) +
    stat_halfeye(
      alpha = 0.7,
      point_interval = median_qi,
      .width = c(.5, .8, .95)
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(
      y = "Predictor",
      x = "Posterior estimate (log scale, negative binomial)",
      fill = NULL
    ) +
    scale_fill_manual(
      values = c(
        "Main effect" = "#3A2767",
        "Interaction" = "#B19CD9"
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom"
    )
  
  
  # Need the centered adult-abundance values to set a sensible range
  dt <- read.csv(dt_file)
  dt$n_adu_log_c <- as.numeric(scale(dt$n_adu_log, scale = FALSE))
  
  adu_grid <- seq(
    min(dt$n_adu_log_c, na.rm = TRUE),
    max(dt$n_adu_log_c, na.rm = TRUE),
    length.out = 60
  )
  
  # For each interacting predictor X, the conditional slope at adult abundance a is:
  # b_X + b_{n_adu_log_c:X} * a
  build_slope <- function(main_col, int_col, label) {
    map_dfr(adu_grid, function(a) {
      slope <- draws[[main_col]] + draws[[int_col]] * a
      
      tibble(
        predictor   = label,
        n_adu_log_c = a,
        slope       = slope
      )
    })
  }
  
  slopes_long <- bind_rows(
    build_slope(
      "b_mean_sett_log_c",
      "b_n_adu_log_c:mean_sett_log_c",
      "average settlers"
    ),
    build_slope(
      "b_R_log10_c",
      "b_n_adu_log_c:R_log10_c",
      "surface rugosity"
    ),
    build_slope(
      "b_Dc",
      "b_n_adu_log_c:Dc",
      "fractal dimension"
    )
  )
  
  # Summarise to median + 50/80/95% intervals
  slopes_summary <- slopes_long %>%
    group_by(predictor, n_adu_log_c) %>%
    median_qi(slope, .width = c(.5, .8, .95)) %>%
    ungroup() %>%
    mutate(
      predictor = factor(
        predictor,
        levels = c(
          "average settlers",
          "surface rugosity",
          "fractal dimension"
        )
      )
    )
  
  cs <- ggplot(
    slopes_summary,
    aes(
      x = n_adu_log_c,
      y = slope,
      ymin = .lower,
      ymax = .upper
    )
  ) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_lineribbon(aes(fill = factor(.width)), alpha = 0.5) +
    geom_line(linewidth = 0.7, colour = "#3A2767") +
    geom_rug(
      data = dt,
      aes(x = n_adu_log_c),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.5
    )+
    facet_wrap(~ predictor, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "0.5"  = "#3A2767",
        "0.8"  = "#7A5BA8",
        "0.95" = "#B19CD9"
      ),
      name = "Credible interval"
    ) +
    labs(
      x = "adult abundance (log, centered)",
      y = "Conditional slopes\nfor Juvenile abundance"
    ) +
    guides(
      fill = guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.margin = margin(t = 4, r = 0, b = 0, l = 0),
      plot.margin = margin(t = 12, r = 8, b = 12, l = 8)
    )
  
  # Free the bottom panel from patchwork alignment with the wide y-axis labels above
  p <- e / patchwork::free(cs) +
    plot_layout(
      heights = c(1, 1.35)
    ) +
    plot_annotation(
      tag_levels = "a"
    ) &
    theme(
      plot.tag = element_text(face = "bold", size = 18),
      plot.tag.position = c(0, 1)
    )
  
  ggsave(
    out_file,
    p,
    width = 8,
    height = 7
  )
  
  out_file_png <- sub("\\.pdf$", ".png", out_file)
  ggsave(
    out_file_png,
    p,
    width = 8,
    height = 7
  )


for(i in seq_len(nrow(scales))) {
  make_interaction_plot(
    mod_file = scales$mod[i],
    dt_file  = scales$dt[i],
    outfile  = scales$out[i]
  )
}

  
####get probabilities ####
  
  interaction_probs <- draws_long %>%
    filter(type == "Interaction") %>%
    group_by(predictor) %>%
    summarise(
      posterior_mean = mean(estimate, na.rm = TRUE),
      median = median(estimate, na.rm = TRUE),
      l95 = quantile(estimate, 0.025, na.rm = TRUE),
      u95 = quantile(estimate, 0.975, na.rm = TRUE),
      p_less_than_0 = mean(estimate < 0, na.rm = TRUE),
      p_greater_than_0 = mean(estimate > 0, na.rm = TRUE),
      probability_of_direction = max(p_less_than_0, p_greater_than_0),
      .groups = "drop"
    )
  
  interaction_probs
  
  
 #### get change in juvs for change in adult ####
  
  beta <- 0.387
  beta_ci <- c(0.31, 0.46)
  
  # Effect of changing from A_old to A_new
  effect_log1p <- function(beta, A_old, A_new) {
    delta_x <- log1p(A_new) - log1p(A_old)
    (exp(beta * delta_x) - 1) * 100
  }
  
  # Doubling effect from A to 2A
  doubling_effect_log1p <- function(beta, A) {
    effect_log1p(beta, A_old = A, A_new = 2 * A)
  }
  
  # Adding one adult
  adding_one_effect_log1p <- function(beta, A) {
    effect_log1p(beta, A_old = A, A_new = A + 1)
  }
  
  # Examples
  A_values <- c(0, 1, 5, 10, 25, 50, 100)
  
  data.frame(
    adults_start = A_values,
    adults_doubled = 2 * A_values,
    doubling_percent_change = doubling_effect_log1p(beta, A_values),
    adults_plus_one = A_values + 1,
    adding_one_percent_change = adding_one_effect_log1p(beta, A_values)
  )
  
  
  