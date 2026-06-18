library(brms)
library(tidyverse)
library(tidybayes)
library(posterior)
library(loo)
library(ggplot2)
library(patchwork)
library(rstan)

# path to your saved models
setwd("~/Documents/LIRSjuvenile_spatial")
model_dir <- "output/models"

# list all models
model_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
model_files
model19_files <- list.files(model_dir, pattern = "\\c19.rds$", full.names = TRUE)
model50_19_files <- model19_files[grepl("dt5", model19_files)]
model1_19_files <- model19_files[grepl("dt1", model19_files)]
model2_19_files <- model19_files[grepl("dt2", model19_files)]
model4_19_files <- model19_files[grepl("dt4", model19_files)]
model23_files <- list.files(model_dir, pattern = "\\c23.rds$", full.names = TRUE)
model50_23_files <- model23_files[grepl("dt5", model23_files)]
model1_23_files <- model23_files[grepl("dt1", model23_files)]
model2_23_files <- model23_files[grepl("dt2", model23_files)]
model4_23_files <- model23_files[grepl("dt4", model23_files)]
modelfull_files <- list.files(model_dir, pattern = "\\_c.rds$", full.names = TRUE)

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


# full models
model5_files <- c("output/models/brm_nb_dt5_full_c.rds","output/models/brm_nb_dt5_R_c.rds", "output/models/brm_nb_dt5_D_c.rds", "output/models/brm_nb_dt5_a_c.rds", "output/models/brm_nb_dt5_s_c.rds")
model1_files <- c("output/models/brm_nb_dt1_full_c.rds","output/models/brm_nb_dt1_R_c.rds", "output/models/brm_nb_dt1_D_c.rds", "output/models/brm_nb_dt1_a_c.rds", "output/models/brm_nb_dt1_s_c.rds")
model2_files <- c("output/models/brm_nb_dt2_full_c.rds","output/models/brm_nb_dt2_R_c.rds", "output/models/brm_nb_dt2_D_c.rds", "output/models/brm_nb_dt2_a_c.rds", "output/models/brm_nb_dt2_s_c.rds")
model4_files <- c("output/models/brm_nb_dt4_full_c.rds","output/models/brm_nb_dt4_R_c.rds", "output/models/brm_nb_dt4_D_c.rds", "output/models/brm_nb_dt4_a_c.rds", "output/models/brm_nb_dt4_s_c.rds")

model_files <- c(model1_files, model2_files, model4_files)

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
#write_csv(model_summary, "output/dt1_cVSs_split_model_comparison_summary.csv")
#write_csv(effects_table, "output/dt1_cVSs_split_model_effects_summary.csv")

effects_clean <- effects_table %>%
  filter(term != "b_Intercept")%>%
  mutate(
    scale = str_sub(model,start=10 , end=10),
    # clean names
    term = str_remove(term, "^b_"),
    term = str_replace_all(term, "_", " "),
    scale = str_replace_all(scale, "1", "scale = 1mx1m"),
    scale = str_replace_all(scale, "2", "scale = 2mx2m"),
    scale = str_replace_all(scale, "4", "scale = 4mx4m"),
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
    term = str_replace_all(term, "Dc", "D"),
    term = str_replace_all(term, "R log10 c","R log10"),
    term = str_replace_all(term, "mean sett log c", "settlers (log)"),
    term = str_replace_all(term, "n adu log c", "adults (log)"),
    

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
  
  # split random effects by group
  if (effect_type == "random") {
    p <- p + facet_wrap(~group, scales = "free_y")
  }
  
  return(p)
}

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
  "figures/subsites_model_effects_caterpillar.png",
  final_plot,
  width = 12,
  height = 30,   # tall because lots of rows
  dpi = 300
)




eff<- effects_clean %>% 
  filter(type == "fixed") %>%
  ggplot(aes(x = mean, y = var, col = var)) +
  
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
  labs(x = NULL, y = "Effect Sizes") +
  facet_grid(term~scale,
             scales = "free",
             switch = "y")+
  scale_colour_manual(
    name = "Fitted model",
    values = c(
      "full model" = "black",
      "alone" = "darkorchid2"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 14),
    strip.text.y.left = element_text(angle = 0, size = 12),
    axis.title = element_text(size = 14),
      strip.placement = "outside",
      strip.background = element_blank()
  )

ggsave(
  "figures/model_effects_caterpillar_color.png",
  eff,
  width = 10,
  height = 4,   # tall because lots of rows
  dpi = 300
)



## effect figure #####

mod <- readRDS("output/models/brm_nb_dt1_full_c.rds")
dt<- read.csv("output/dt1.csv")

dt$site <- as.factor(dt$site)
dt$year <- as.factor(dt$year)
dt$n_adu_log_c     <- scale(dt$n_adu_log, scale = FALSE)
dt$mean_sett_log_c <- scale(dt$mean_sett_log, scale = FALSE)
dt$R_log10_c       <- scale(dt$R_log10, scale = FALSE)
dt$Dc            <- scale(dt$D, scale = FALSE)
dt$n_adu_log_s     <- scale(dt$n_adu_log)
dt$mean_sett_log_s <- scale(dt$mean_sett_log)
dt$R_log10_s       <- scale(dt$R_log10)
dt$Ds            <- scale(dt$D)


# adults
ma <-readRDS("output/models/brm_nb_dt1_a_c.rds")
pred_adu <- data.frame(
  n_adu_log_c = seq(min(dt$n_adu_log_c),max(dt$n_adu_log_c),length.out = 500),
  R_log10_c = 0,
  Dc = 0,
  mean_sett_log_c = 0)

fit_adu <- fitted(
  mod,
  newdata = pred_adu,
  re_formula = NA,
  summary = TRUE
)
fit_adu2 <- fitted(
  ma,
  newdata = pred_adu,
  re_formula = NA,
  summary = TRUE
)

pred_adu$estimate <- fit_adu[, "Estimate"]
pred_adu$lower    <- fit_adu[, "Q2.5"]
pred_adu$upper    <- fit_adu[, "Q97.5"]

pred_adu$estimate2 <- fit_adu2[, "Estimate"]
pred_adu$lower2    <- fit_adu2[, "Q2.5"]
pred_adu$upper2    <- fit_adu2[, "Q97.5"]

a <- ggplot() +
  # observed data
  geom_jitter(
    data = dt,
    aes(x = n_adu_log_c, y = n_juv),
    alpha = 0.08,
    size = 1.8,
    col = "darkorchid4"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_adu,
    aes(
      x = n_adu_log_c,
      ymin = lower,
      ymax = upper
    ),
    #fill = "darkorchid1",
    alpha = 0.2
  ) +
  
  # fitted model
  geom_line(
    data = pred_adu,
    aes(x = n_adu_log_c, y = estimate),
    linewidth = 1.1,
    alpha = 0.7
    #col = "darkorchid1"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_adu,
    aes(
      x = n_adu_log_c,
      ymin = lower2,
      ymax = upper2
    ),
    fill = "orchid",
    alpha = 0.1
  ) +
  theme_classic(base_size = 14) +
  scale_colour_manual(
    name = "Fitted model",
    values = c(
      "full model" = "black",
      "one predictor model" = "orchid"
    )
  ) +
  coord_cartesian(ylim = c(0, 120))+
  # fitted model
  geom_line(
    data = pred_adu,
    aes(x = n_adu_log_c, y = estimate2),
    linewidth = 1.1,
    alpha = 0.7,
    col = "orchid"
  ) +
  
  
  
  labs(
    x = "Adult abundance (log-centered)",
    y = "Juvenile abundance"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(colour = "black"),
    panel.border = element_blank(),
    legend.position = "none"
  )


# R
mr <-readRDS("output/models/brm_nb_dt1_R_c.rds")
pred_R <- data.frame(
  n_adu_log_c = 0,
  R_log10_c = seq(min(dt$R_log10_c),max(dt$R_log10_c),length.out = 500),
  Dc = 0,
  mean_sett_log_c = 0)

fit_r <- fitted(
  mod,
  newdata = pred_R,
  re_formula = NA,
  summary = TRUE
)
fit_r2 <- fitted(
  mr,
  newdata = pred_R,
  re_formula = NA,
  summary = TRUE
)

pred_R$estimate <- fit_r[, "Estimate"]
pred_R$lower    <- fit_r[, "Q2.5"]
pred_R$upper    <- fit_r[, "Q97.5"]

pred_R$estimate2 <- fit_r2[, "Estimate"]
pred_R$lower2    <- fit_r2[, "Q2.5"]
pred_R$upper2    <- fit_r2[, "Q97.5"]

r <- ggplot() +
  # observed data
  geom_point(
    data = dt,
    aes(x = R_log10_c, y = n_juv),
    alpha = 0.08,
    size = 1.8,
    col = "darkorchid4"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_R,
    aes(
      x = R_log10_c,
      ymin = lower,
      ymax = upper
    ),
    #fill = "darkorchid1",
    alpha = 0.2
  ) +
  
  # fitted model
  geom_line(
    data = pred_R,
    aes(x = R_log10_c, y = estimate),
    linewidth = 1.1,
    alpha = 0.7
    #col = "darkorchid1"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_R,
    aes(
      x = R_log10_c,
      ymin = lower2,
      ymax = upper2
    ),
    fill = "orchid",
    alpha = 0.1
  ) +
  
  # fitted model
  geom_line(
    data = pred_R,
    aes(x = R_log10_c, y = estimate2),
    linewidth = 1.1,
    alpha = 0.7,
    col = "orchid"
  ) +
  theme_classic(base_size = 14) +
  scale_colour_manual(
    name = "Fitted model",
    values = c(
      "full model" = "black",
      "one predictor model" = "orchid"
    )
  ) +
  coord_cartesian(ylim = c(0, 120))+
  labs(
    x = "Surface rugosity (log-centered)",
    y = "Juvenile abundance"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(colour = "black"),
    panel.border = element_blank(),
    legend.position = "none"
  )



# D
md <-readRDS("output/models/brm_nb_dt1_D_c.rds")
pred_D <- data.frame(
  n_adu_log_c = 0,
  R_log10_c = 0,
  Dc = seq(min(dt$Dc),max(dt$Dc),length.out = 500),
  mean_sett_log_c = 0)

fit_D <- fitted(
  mod,
  newdata = pred_D,
  re_formula = NA,
  summary = TRUE
)
fit_D2 <- fitted(
  md,
  newdata = pred_D,
  re_formula = NA,
  summary = TRUE
)

pred_D$estimate <- fit_D[, "Estimate"]
pred_D$lower    <- fit_D[, "Q2.5"]
pred_D$upper    <- fit_D[, "Q97.5"]

pred_D$estimate2 <- fit_D2[, "Estimate"]
pred_D$lower2    <- fit_D2[, "Q2.5"]
pred_D$upper2    <- fit_D2[, "Q97.5"]

d <- ggplot() +
  # observed data
  geom_point(
    data = dt,
    aes(x = Dc, y = n_juv),
    alpha = 0.08,
    size = 1.8,
    col = "darkorchid4"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_D,
    aes(
      x = Dc,
      ymin = lower,
      ymax = upper
    ),
    #fill = "darkorchid1",
    alpha = 0.2
  ) +
  
  # fitted model
  geom_line(
    data = pred_D,
    aes(x = Dc, y = estimate),
    linewidth = 1.1,
    alpha = 0.7
    #col = "darkorchid1"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_D,
    aes(
      x = Dc,
      ymin = lower2,
      ymax = upper2
    ),
    fill = "orchid",
    alpha = 0.1
  ) +
  
  # fitted model
  geom_line(
    data = pred_D,
    aes(x = Dc, y = estimate2),
    linewidth = 1.1,
    alpha = 0.7,
    col = "orchid"
  ) +
  theme_classic(base_size = 14) +
  scale_colour_manual(
    name = "Fitted model",
    values = c(
      "full model" = "black",
      "one predictor model" = "orchid"
    )
  ) +
  coord_cartesian(ylim = c(0, 120))+
  labs(
    x = "Fractal dimension (centered)",
    y = "Juvenile abundance"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(colour = "black"),
    panel.border = element_blank(),
    legend.position = "none"
  )

  
# S
ms <-readRDS("output/models/brm_nb_dt1_s_c.rds")
pred_s19 <- data.frame(
  n_adu_log_c = 0,
  R_log10_c = 0,
  Dc = 0,
  mean_sett_log_c = seq(min(dt$mean_sett_log_c),max(dt$mean_sett_log_c ),length.out = 500),
  year = as.factor("2019")
  )

fit_s <- fitted(
  mod,
  newdata = pred_s,
  re_formula = NA,
  summary = TRUE
)
fit_s2 <- fitted(
  ms,
  newdata = pred_s,
  re_formula = NA,
  summary = TRUE
)

fit_s19 <- fitted(
  ms,
  newdata = pred_s,
  re_formula = ~(1|year),
  summary = TRUE
)

pred_s23 <- data.frame(
  n_adu_log_c = 0,
  R_log10_c = 0,
  Dc = 0,
  mean_sett_log_c = seq(min(dt$mean_sett_log_c),max(dt$mean_sett_log_c ),length.out = 500),
  year = as.factor("2023")
)

fit_s23 <- fitted(
  ms,
  newdata = pred_s23,
  re_formula = ~(1|year),
  summary = TRUE
)

pred_s$estimate <- fit_s[, "Estimate"]
pred_s$lower    <- fit_s[, "Q2.5"]
pred_s$upper    <- fit_s[, "Q97.5"]

pred_s$estimate19 <- fit_s19[, "Estimate"]
pred_s$lower19    <- fit_s19[, "Q2.5"]
pred_s$upper19    <- fit_s19[, "Q97.5"]

pred_s$estimate2 <- fit_s2[, "Estimate"]
pred_s$lower2    <- fit_s2[, "Q2.5"]
pred_s$upper2    <- fit_s2[, "Q97.5"]

pred_s$estimate23 <- fit_s23[, "Estimate"]
pred_s$lower23   <- fit_s23[, "Q2.5"]
pred_s$upper23    <- fit_s23[, "Q97.5"]


s2 <- ggplot() +
  # observed data
  geom_jitter(
    data = dt,
    aes(x = mean_sett_log_c, y = n_juv,   col = year),
    alpha = 0.08,
    size = 1.8
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_s,
    aes(
      x = mean_sett_log_c,
      ymin = lower,
      ymax = upper
    ),
    #fill = "darkorchid1",
    alpha = 0.2
  ) +
  
  # fitted model
  geom_line(
    data = pred_s,
    aes(x = mean_sett_log_c, y = estimate),
    linewidth = 1.1,
    alpha = 0.7
    #col = "darkorchid1"
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_s,
    aes(
      x = mean_sett_log_c,
      ymin = lower2,
      ymax = upper2
    ),
    fill = "orchid",
    alpha = 0.1
  ) +
  
  # fitted model
  geom_line(
    data = pred_s,
    aes(x = mean_sett_log_c, y = estimate2),
    linewidth = 1.1,
    alpha = 0.7,
    col = "orchid"
  ) +
  
  geom_line(
    data = pred_s,
    aes(x = mean_sett_log_c, y = estimate19),
    linewidth = 1.1,
    alpha = 0.7,
    col = "blue"
  ) +
  
  geom_line(
    data = pred_s,
    aes(x = mean_sett_log_c, y = estimate23),
    linewidth = 1.1,
    alpha = 0.7,
    col = "green"
  ) +
  
  
  labs(
    x = "mean settlers (log-centered)",
    y = "Juvenile abundance"
  ) +
  
  theme_classic(base_size = 14) +
  # scale_colour_manual(
  #   name = "Fitted model",
  #   values = c(
  #     "full model" = "black",
  #     "one predictor model" = "orchid"
  #   )
  # ) +
  coord_cartesian(ylim = c(0, 120))+
  theme(
    axis.text = element_text(colour = "black"),
    panel.border = element_blank(),
    legend.position = "none"
  )

((a + d) /
    (r + s)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



### R vs D


pred_grid <- expand.grid(
  R_log10_c = seq(
    min(dt$R_log10_c),
    max(dt$R_log10_c),
    length.out = 500
  ),
  Dc = seq(
    min(dt$Dc),
    max(dt$Dc),
    length.out = 500
  )
)

pred_grid$n_adu_log_c    <- 0
pred_grid$mean_sett_log_c <- 0

fit <- fitted(
  mod,
  newdata = pred_grid,
  re_formula = NA,
  summary = TRUE
)

pred_grid$estimate <- fit[, "Estimate"]
pred_grid$lower    <- fit[, "Q2.5"]
pred_grid$upper    <- fit[, "Q97.5"]

ggplot(pred_grid,
       aes(R_log10_c, Dc, fill = estimate)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Predicted\nn_juv") +
  coord_equal() +
  theme_bw()


library(sf)

# observed points
obs <- st_as_sf(
  dt,
  coords = c("R_log10_c", "Dc"),
  crs = NA
)

# convex hull
hull <- obs |>
  st_union() |>
  st_convex_hull()

# prediction points
pred_pts <- st_as_sf(
  pred_grid,
  coords = c("R_log10_c", "Dc"),
  crs = NA
)

# keep only points inside hull
inside <- st_within(pred_pts, hull, sparse = FALSE)[,1]

pred_grid$estimate_clip <- pred_grid$estimate
pred_grid$estimate_clip[!inside] <- NA

ggplot(pred_grid,
       aes(R_log10_c, Dc, fill = estimate_clip, alpha = 0.3)) +
  geom_raster() +
  geom_point(
    data = dt,
    aes(x = R_log10_c, y = Dc, col = n_juv_log),
    inherit.aes = FALSE,
    pch = 16,
    size = 3,
    alpha = 0.8
  ) +
  scale_fill_viridis_c(
    name = "Predicted\nn_juv",
    na.value = NA
  ) +
  theme_minimal()
