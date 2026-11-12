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

model5_files <- c("output/models/brm_nb_dt5_full_c_aint.rds","output/models/brm_nb_dt5_full_c_int.rds", "output/models/brm_nb_dt5_full_c.rds")
model1_files <- c("output/models/brm_nb_dt1_full_c_aint.rds","output/models/brm_nb_dt1_full_c_int.rds", "output/models/brm_nb_dt1_full_c.rds")
model2_files <- c("output/models/brm_nb_dt2_full_c_aint.rds","output/models/brm_nb_dt2_full_c_int.rds", "output/models/brm_nb_dt2_full_c.rds")
model4_files <- c("output/models/brm_nb_dt4_full_c_aint.rds","output/models/brm_nb_dt4_full_c_int.rds", "output/models/brm_nb_dt4_full_c.rds")


model_files <- c(model1_files, model2_files, model5_files)


model_files <- c("output/models/brm_nb_dt5_full_c_aint_ry.rds", "output/models/brm_nb_dt5_full_c_aint.rds",
                 "output/models/brm_nb_dt1_full_c_aint_ry.rds", "output/models/brm_nb_dt1_full_c_aint.rds",
                 "output/models/brm_nb_dt2_full_c_aint_ry.rds", "output/models/brm_nb_dt2_full_c_aint.rds",
                 "output/models/brm_nb_dt4_full_c_aint_ry.rds","output/models/brm_nb_dt4_full_c_aint.rds" )
                 "output/models/brm_nb_dt1_full_c.rds", "output/models/brm_nb_dt1_full_c_int.rds",
                 "output/models/brm_nb_dt2_full_c.rds", "output/models/brm_nb_dt2_full_c_int.rds",
                 "output/models/brm_nb_dt4_full_c.rds", "output/models/brm_nb_dt4_full_c_int.rds")
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

title_df <- tibble::tribble(
  ~model, ~title,
  "brm_nb_dt5_full_c", "50x50cm scale, year as random effect",
  "brm_nb_dt5_full_c_int", "50x50cm scale, year as interacting effect",
  "brm_nb_dt5_full_c_aint", "50x50cm scale, adu abundance as interacting effect",
  "brm_nb_dt5_full_c_aint_ry", "50x50cm scale, adu abu interaction, year random",
  "brm_nb_dt1_full_c", "1x1m scale, year as random effect",
  "brm_nb_dt1_full_c_int", "1x1m scale, year as interacting effect",
  "brm_nb_dt1_full_c_aint", "1x1m scale, adu abundance as interacting effect",
  "brm_nb_dt1_full_c_aint_ry", "1x1m scale, adu abu interaction, year random",
  "brm_nb_dt2_full_c", "2x2m scale, year as random effect",
  "brm_nb_dt2_full_c_int", "2x2m scale, year as interacting effect",
  "brm_nb_dt2_full_c_aint", "2x2m scale, adu abundance as interacting effect",
  "brm_nb_dt2_full_c_aint_ry", "2x2m scale, adu abu interaction, year random",
  "brm_nb_dt4_full_c", "4x4m scale, year as random effect",
  "brm_nb_dt4_full_c_int", "4x4m scale, year as interacting effect",
  "brm_nb_dt4_full_c_aint", "4x4m scale, adu abundance as interacting effect",
  "brm_nb_dt4_full_c_aint_ry", "4x4m scale, adu abu interaction, year random"
)

plot_rows <- lapply(models, function(m) {
  
  df_m <- effects_clean %>%
    filter(model == m)
  
  plot_title <- title_df %>%
    filter(model == m) %>%
    pull(title)
  
  p_fixed <- plot_effects(df_m, "fixed") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(
        size = 10,
        face = "bold")
    )
  
  p_random <- plot_effects(df_m, "random")
  
  p_fixed + p_random +
    plot_layout(widths = c(1, 1.2))
})

final_plot <- wrap_plots(plot_rows, ncol = 2)
final_plot

ggsave("figures/models_effects_caterpillar.png",
  final_plot,width = 12,height = 10,dpi = 300
)




## to run with alone and full models

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
### interaction model #####

make_interaction_plot <- function(mod_file, dt_file, outfile){
  
  mod <- readRDS(mod_file)
  dt  <- read.csv(dt_file)
  dt$site <- factor(dt$site)
  dt$year <- factor(dt$year)
  dt$n_adu_log_c     <- scale(dt$n_adu_log, scale = FALSE)
  dt$mean_sett_log_c <- scale(dt$mean_sett_log, scale = FALSE)
  dt$R_log10_c       <- scale(dt$R_log10, scale = FALSE)
  dt$Dc              <- scale(dt$D, scale = FALSE)
  
  make_panel <- function(var, xlabel){
    
    pred_grid <- expand.grid(year = levels(dt$year),
      value = seq( min(dt[[var]], na.rm = TRUE),max(dt[[var]], na.rm = TRUE),length.out = 500))
    
    for (p in preds){pred_grid[[p]] <- 0}
    pred_grid[[var]] <- pred_grid$value
    pred_grid$value <- NULL
    
    fit <- fitted(mod,newdata = pred_grid, re_formula = NA,summary = TRUE)
    
    pred_grid$estimate <- fit[, "Estimate"]
    pred_grid$lower    <- fit[, "Q2.5"]
    pred_grid$upper    <- fit[, "Q97.5"]
    
    ggplot() +
      geom_ribbon(data = pred_grid,
        aes_string(x = var, fill = "year", ymin = "lower", ymax = "upper"),
        alpha = 0.2) +
      geom_jitter(data = dt,
        aes_string(x = var,y = "n_juv",colour = "year"),
        alpha = 0.12,size = 1.8) +
      geom_line(data = pred_grid,
        aes_string(x = var,y = "estimate",colour = "year"),
        linewidth = 1.1,alpha = 0.7) +
      theme_classic(base_size = 14) +
      scale_color_manual(values = c("#B19CD9", "#3A2767")) +
      scale_fill_manual(values = c("#B19CD9", "#3A2767")) +
      labs(x = xlabel, y = "Juvenile abundance") +
      theme(
        axis.text = element_text(colour = "black"),
        panel.border = element_blank())
  }
  
  a_int <- make_panel( "n_adu_log_c","Adult abundance (log, centered)")
  d_int <- make_panel("Dc", "Fractal dimension")
  r_int <- make_panel("R_log10_c","Surface rugosity (log10, centered)")
  s_int <- make_panel("mean_sett_log_c","Average settlers (log, centered)")
  
  p <- ((a_int + d_int) /(r_int + s_int))/(e+ d)+
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave( outfile,p, width = 10,height = 8)
  
  return(p)
}


scales <- data.frame(
  mod = c( "output/models/brm_nb_dt1_full_c_int.rds","output/models/brm_nb_dt2_full_c_int.rds",
    "output/models/brm_nb_dt4_full_c_int.rds","output/models/brm_nb_dt5_full_c_int.rds"),
  dt = c("output/dt1.csv","output/dt2.csv","output/dt4.csv","output/dt50.csv"),
  out = c("figures/interaction_panels_dt1.pdf","figures/interaction_panels_dt2.pdf",
    "figures/interaction_panels_dt4.pdf","figures/interaction_panels_dt50.pdf"))

preds <- c("n_adu_log_c","mean_sett_log_c","R_log10_c","Dc")

for(i in seq_len(nrow(scales))) {
  make_interaction_plot(
    mod_file = scales$mod[i],
    dt_file  = scales$dt[i],
    outfile  = scales$out[i]
  )
}

mod_list <- c("brm_nb_dt1_full_c_int","brm_nb_dt2_full_c_int",
              "brm_nb_dt4_full_c_int","brm_nb_dt5_full_c_int")


for(m in mod_list) {
df_m <- effects_clean %>%
  filter(model == m) %>%
  mutate(
    term_label = recode(
      term,
      "adults (log)" = "Adult abundance",
      "settlers (log)" = "Average settlers",
      "R log10" = "Surface rugosity",
      "D" = "Fractal dimension",
      "year2023:adults (log)" = "Year × Adult abundance",
      "year2023:settlers (log)" = "Year × Average settlers",
      "year2023:R log10" = "Year × Surface rugosity",
      "year2023:D" = "Year × Fractal dimension",
      "year2023" = "Year (2023)"
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
  geom_errorbarh(aes(xmin = l95, xmax = u95),
                 height = 0, size = 0.4, alpha = 0.7) +
  geom_errorbarh(aes(xmin = l80, xmax = u80),
                 height = 0, size = 1.2) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  scale_x_continuous(limits = x_limits) +
  labs(x = "Model effect sizes", y = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
  )

mod <- readRDS(paste0("output/models/", m, ".rds"))
draws <- as_draws_df(mod)
draws_long <- tibble(
  `Adult abundance_2019` = draws$b_n_adu_log_c,
  `Adult abundance_2023` = draws$b_n_adu_log_c + draws$`b_year2023:n_adu_log_c`,
  `Average settlers_2019` = draws$b_mean_sett_log_c,
  `Average settlers_2023` = draws$b_mean_sett_log_c + draws$`b_year2023:mean_sett_log_c`,
  `Surface rugosity_2019` = draws$b_R_log10_c,
  `Surface rugosity_2023` = draws$b_R_log10_c + draws$`b_year2023:R_log10_c`,
  `Fractal dimension_2019` = draws$b_Dc,
  `Fractal dimension_2023` = draws$b_Dc + draws$`b_year2023:Dc`) %>%
  pivot_longer(everything(),
    names_to = c("predictor", "year"), names_sep = "_", values_to = "estimate")

draws_long$predictor <- factor(
  draws_long$predictor,
  levels = c("Adult abundance", "Average settlers", "Surface rugosity", "Fractal dimension"))

d <- ggplot(draws_long, aes(x = estimate,y = predictor,fill = year)) +
  stat_halfeye(alpha = 0.6, point_interval = median_qi,.width = c(.5, .8, .95),position = position_dodge(width = 0.4)) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(y="Predictor", x = "Posterior draws estimates")+
  scale_fill_manual(values = c("2019" = "#B19CD9","2023" = "#3A2767")) +
  theme_classic(base_size = 14)

p <- e + d 
ggsave(paste0("figures/",m,"_effSizes_and_postDraws.pdf"),p, width = 8,height = 3.5)

}





scales <- data.frame(
  mod = c( "output/models/brm_nb_dt1_full_s_int.rds","output/models/brm_nb_dt2_full_s_int.rds",
           "output/models/brm_nb_dt4_full_s_int.rds","output/models/brm_nb_dt5_full_s_int.rds"),
  dt = c("output/dt1.csv","output/dt2.csv","output/dt4.csv","output/dt50.csv"),
  out = c("figures/interaction_panels_dt1_s.pdf","figures/interaction_panels_dt2_s.pdf",
          "figures/interaction_panels_dt4_s.pdf","figures/interaction_panels_dt50_s.pdf"))

preds <- c("n_adu_log_s","mean_sett_log_s","R_log10_s","Ds")

for(i in seq_len(nrow(scales))) {
  make_interaction_plot(
    mod_file = scales$mod[i],
    dt_file  = scales$dt[i],
    outfile  = scales$out[i]
  )
}


