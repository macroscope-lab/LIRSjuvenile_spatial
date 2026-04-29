# libraries

library(brms)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(patchwork)
library(posterior)

# read data ####

paths <- c("output/dt50.csv","output/dt1.csv","output/dt2.csv",
           "output/dt4.csv")


## functions for plots #####

coef_plot <- function(fit) {
  
  draws <- as_draws_df(fit) %>% as_tibble()
  
  # fixed effects
  fix_df <- draws %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      everything(),
      names_to = "term",
      values_to = "value"
    ) %>%
    mutate(
      term = str_remove(term, "^b_"),
      group = "Fixed"
    )
  
  # random effects
  re_df_site <- draws %>%
    select(matches("^r_site\\[.*,Intercept\\]$")) %>%
    pivot_longer(
      everything(),
      names_to = "term",
      values_to = "value"
    ) %>%
    mutate(
      term = str_extract(term, "(?<=\\[).*(?=,)"),
      term = paste0("site: ", term),
      group = "Random"
    )
  
  re_df_year<- draws %>%
    select(matches("^r_year\\[.*,Intercept\\]$")) %>%
    pivot_longer(
      everything(),
      names_to = "term",
      values_to = "value"
    ) %>%
    mutate(
      term = str_extract(term, "(?<=\\[).*(?=,)"),
      term = paste0("year: ", term),
      group = "Random"
    )
  
  all_df <- bind_rows(fix_df, re_df_site, re_df_year)
  
  sum_df <- all_df %>%
    group_by(group, term) %>%
    median_qi(value, .width = c(.90, .95)) %>%
    ungroup()
  
  ggplot(sum_df,
         aes(x = value,
             y = fct_reorder(term, value))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_linerange(
      data = filter(sum_df, .width == 0.95),
      aes(xmin = .lower, xmax = .upper),
      linewidth = 1.4
    ) +
    geom_linerange(
      data = filter(sum_df, .width == 0.90),
      aes(xmin = .lower, xmax = .upper),
      linewidth = 0.6
    ) +
    geom_point(size = 2) +
    facet_wrap(~group, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    
    labs(x = "Estimate", y = NULL)
}

obs_pred_plot <- function(fit, dt) {
  dt$pred <- fitted(fit)[, "Estimate"]
  ggplot(dt, aes(x = pred, y = n_juv)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw(base_size = 12) +
    labs(
      x = "Predicted",
      y = "Observed",
      title = "Observed vs Predicted"
    )
}


for (path in paths) {
  dt <- read.csv(path)
  dt$site <- as.factor(dt$site)
  dt$year <- as.factor(dt$year)
  
  fitnb_full <- brm(
    n_juv ~ n_adu_log + mean_sett_log + R_log10 + Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_full,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_full.rds"))
  
  fitnb_inter <- brm(
    n_juv ~ n_adu_log*(mean_sett_log + R_log10 + Ds) + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_inter,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_full_inter.rds"))
  
  fitnb_nosett <- brm(
    n_juv ~ n_adu_log + R_log10 + Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_nosett,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_nosett.rds"))
  
  fitnb_Ra <- brm(
    n_juv ~ n_adu_log + R_log10 + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_Ra,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_Ra.rds"))
  
  fitnb_Ra_int <- brm(
    n_juv ~ n_adu_log * R_log10 + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_Ra_int,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_Ra_int.rds"))
  
  fitnb_Da <- brm(
    n_juv ~ n_adu_log + Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_Da,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_Da.rds"))
  
  fitnb_Da_int <- brm(
    n_juv ~ n_adu_log * Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_Da_int,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_Da_int.rds"))
  
  fitnb_DR <- brm(
    n_juv ~ R_log10 + Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.90, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_DR,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_DR.rds"))
  
  fitnb_DR_int <- brm(
    n_juv ~ R_log10 + Ds + (1|year) + (1 | site),
    data = dt,
    iter = 10000, warmup = 2000, thin = 5,
    control = list(adapt_delta = 0.98, max_treedepth = 15),
    family = negbinomial())
  saveRDS(fitnb_DR_int,paste0("output/models/brm_nb_",str_sub(path, 8, 10),"_DR_int.rds"))
}




  ggsave(
    file.path("figures", paste0("coef_nb_", str_sub(path, 8, 10), ".png")),
    coef_plot(fitnb),
    width = 10, height = 12, dpi = 300
  )
  
  png(file.path("figures", paste0("ppcheck_nb_", str_sub(path, 8, 10), ".png")),
      width = 1800, height = 1200, res = 220)
  print(pp_check(fitnb, ndraws = 100))
  dev.off()
  
  ggsave(
    file.path("figures", paste0("obspred_nb_", str_sub(path, 8, 10), ".png")),
    obs_pred_plot(fitnb, dt),
    width = 8, height = 6, dpi = 300
  )
  
  # fitznb <- brm(
  #   n_juv ~ n_adu_log + mean_sett_log + R_log10 + Ds + (1|year) + (1 | site),
  #   data = dt,
  #   family = zero_inflated_negbinomial(),
  #   prior = c(
  #     prior(normal(0, 5), class = "b"),      # Weakly informative priors for slopes
  #     prior(exponential(1), class = "sd")    # Prior for site-to-site variability
  #   )
  # )
  # saveRDS(fitznb,paste0("output/models/brm_znb_",str_sub(path, 8, 10),".rds"))
  # 
  # ggsave(
  #   file.path("figures", paste0("coef_znb_", str_sub(path, 8, 10), ".png")),
  #   coef_plot(fitznb),
  #   width = 10, height = 12, dpi = 300
  # )
  # 
  # png(file.path("figures", paste0("ppcheck_znb_", str_sub(path, 8, 10), ".png")),
  #     width = 1800, height = 1200, res = 220)
  # print(pp_check(fitznb, ndraws = 100))
  # dev.off()
  # 
  # ggsave(
  #   file.path("figures", paste0("obspred_znb_", str_sub(path, 8, 10), ".png")),
  #   obs_pred_plot(fitznb, dt),
  #   width = 8, height = 6, dpi = 300
  # )
  
  
}






  p <- plot(conditional_effects(fitnb), points = TRUE, ask = FALSE, page = 1)
  
  
  newdata <- data.frame(
    n_adu = seq(min(dt$n_adu), max(dt$n_adu), length.out = 100),
    mean_sett = 0,
    R= 0,
    D = 0,
    site = "horseshoe"  # typical site level
  )
  
  # Get predicted values with uncertainty
  preds <- posterior_epred(fitnb, newdata = newdata, allow_new_levels = TRUE)
  preds2 <- posterior_predict(fitnb, newdata = newdata, allow_new_levels = TRUE)
  

  # Summarize
  pred_mean <- apply(preds, 2, mean)
  pred_lower <- apply(preds, 2, quantile, probs = 0.05)
  pred_upper <- apply(preds, 2, quantile, probs = 0.95)
  
  # Plot
  pred <- ggplot(newdata, aes(x = n_adu, y = pred_mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.2) +
    labs(y = "Predicted juvenile count", x = "Adult count")+
    theme_minimal()
  ggsave(preds, paste0("output/plot_predict_",str_sub(path, 8, 10),".rds"))
  
  # Predicted means
  fitted_vals <- fitted(fitnb)
  
  # Plot predicted vs. observed
  ggplot(data = dt1, aes(y = fitted_vals[, "Estimate"], x = n_juv, col = site)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(y = "Predicted juvenile count", x = "Observed juvenile count") +
    theme_minimal()
  
  ggplot(data = dt4, aes(y = n_juv, x = R, col = site)) +
    geom_point(alpha = 0.4) +
    geom_smooth()+
    geom_smooth(aes(y = n_juv, x = R), col = "black")+
    labs(y = "juvenile count", x = "R") +
    theme_minimal()
  
  ggplot(data = dt4, aes(y = n_juv, x = D, col = site)) +
    geom_point(alpha = 0.4) +
    geom_smooth()+
    geom_smooth(aes(y = n_juv, x = D), col = "black")+
    labs(y = "juvenile count", x = "D") +
    theme_minimal()
  
  
  ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
  site_re_df <- as.data.frame(ranef_vals)
  site_re_df$site <- rownames(ranef_vals)
  site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]
  
  ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
    geom_point() +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(x = "Site", y = "Random Intercept Estimate", 
         title = "Site Random Effects") +
    theme_minimal()
  
  
  fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
  fixef_df$term <- rownames(fixef_df)
  fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]
  
  ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
    geom_point() +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(x = "Fixed Effect", y = "Estimate",
         title = "Fixed Effects (Posterior Means & 95% CI)") +
    theme_minimal()
  
  
  bayesplot::ppc_dens_overlay(y = dt1$n_juv, yrep = posterior_predict(fitnb, draws = 100))
  
}



### 2022 1x1m ####
dt1 <- read.csv("data/dt1_2022.csv")
dt1$site <- as.factor(dt1$site)
dt1$year <- as.factor("2022")


# delete GT square with depth issues
if (unique(dt1$year) == "2022") 
  {#dt1 %>% arrange(desc(H))
  dt1 <- dt1[dt1$H != max(dt1$H),]}
 
# dt1 %>% group_by(site,year)%>%mutate(av_R = mean(R)) %>% select(site, year,av_R) %>% unique() %>% arrange(year,av_R)
# dt1 %>% group_by(site,year)%>%mutate(av_adu = mean(n_adu)) %>% select(site, year,av_adu) %>% unique() %>% arrange(year,av_adu)
# dt1 %>% group_by(site,year)%>%mutate(av_juv = mean(n_juv)) %>% select(site, year,av_juv) %>% unique() %>% arrange(year,av_juv)
dt1 %>% arrange(n_juv) %>% tail(n=10)

ggplot(dt1, aes(n_juv)) +
  geom_histogram(binwidth = 1)+
  theme_bw()+
  ggtitle("cell count per juvenile abundance - 1x1m, 2022")

summary_data <- dt1 %>%
  group_by(site) %>%
  summarise(mean_count = mean(n_juv),
            var_count = var(n_juv))

plot(summary_data$mean_count, summary_data$var_count,
     xlab = "Mean", ylab = "Variance", main = "Mean-Variance across sites - 1x1m, 2022")

poisson_model <- glmer(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s +
                         (1 | site), 
                       data = dt1,
                       family = poisson)

fitznb <- brm(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s +
              (1 | site), 
            data = dt1,
            family = zero_inflated_negbinomial())
fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
                (1 | site), 
              data = dt1,
              family = negbinomial())
fitznb_nat <- brm(n_juv ~ n_adu + mean_sett + R + D +
                (1 | site), 
              data = dt1,
              family = zero_inflated_negbinomial())

fithnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
              (1 | site), 
            data = dt1,
            family = hurdle_negbinomial())
loo(fitnb,fitznb_nat)
plot(conditional_effects(fitznb_nat, effects = "n_adu"))
plot(conditional_effects(fitznb, effects = "n_adu_log"))








fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt1,
             family = negbinomial())




###2022 2x2m ####


fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt2,
             family = negbinomial())



ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
site_re_df <- as.data.frame(ranef_vals)
site_re_df$site <- rownames(ranef_vals)
site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]

ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Site", y = "Random Intercept Estimate", 
       title = "Site Random Effects") +
  theme_minimal()


fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
fixef_df$term <- rownames(fixef_df)
fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]

ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Fixed Effect", y = "Estimate",
       title = "Fixed Effects (Posterior Means & 95% CI)") +
  theme_minimal()


###2022 4x4m ####


fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt4,
             family = negbinomial())



ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
site_re_df <- as.data.frame(ranef_vals)
site_re_df$site <- rownames(ranef_vals)
site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]

ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Site", y = "Random Intercept Estimate", 
       title = "Site Random Effects") +
  theme_minimal()


fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
fixef_df$term <- rownames(fixef_df)
fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]

ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Fixed Effect", y = "Estimate",
       title = "Fixed Effects (Posterior Means & 95% CI)") +
  theme_minimal()


### 2019 1x1m ####
dt1 <- read.csv("data/dt1_2019.csv")
dt1$site <- as.factor(dt1$site)
dt1$year <- as.factor("2019")

fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt1,
             family = negbinomial())



ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
site_re_df <- as.data.frame(ranef_vals)
site_re_df$site <- rownames(ranef_vals)
site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]

ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Site", y = "Random Intercept Estimate", 
       title = "Site Random Effects") +
  theme_minimal()


fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
fixef_df$term <- rownames(fixef_df)
fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]

ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Fixed Effect", y = "Estimate",
       title = "Fixed Effects (Posterior Means & 95% CI)") +
  theme_minimal()




### 2019 2x2m ####
dt2 <- read.csv("data/dt2_2019.csv")
dt2$site <- as.factor(dt2$site)
dt2$year <- as.factor("2019")

fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt2,
             family = negbinomial())



ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
site_re_df <- as.data.frame(ranef_vals)
site_re_df$site <- rownames(ranef_vals)
site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]

ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Site", y = "Random Intercept Estimate", 
       title = "Site Random Effects") +
  theme_minimal()


fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
fixef_df$term <- rownames(fixef_df)
fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]

ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Fixed Effect", y = "Estimate",
       title = "Fixed Effects (Posterior Means & 95% CI)") +
  theme_minimal()


###2019 4x4m ####
dt4 <- read.csv("data/dt4_2019.csv")
dt4$site <- as.factor(dt4$site)
dt4$year <- as.factor("2019")

fitnb <- brm(n_juv ~ n_adu + mean_sett + R + D +
               (1 | site), 
             data = dt4,
             family = negbinomial())



ranef_vals <- ranef(fitnb)$site[, , "Intercept"]
site_re_df <- as.data.frame(ranef_vals)
site_re_df$site <- rownames(ranef_vals)
site_re_df <- site_re_df[, c("site", "Estimate", "Q2.5", "Q97.5")]

ggplot(site_re_df, aes(x = reorder(site, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Site", y = "Random Intercept Estimate", 
       title = "Site Random Effects") +
  theme_minimal()


fixef_df <- as.data.frame(fixef(fitnb))[2:5,]
fixef_df$term <- rownames(fixef_df)
fixef_df <- fixef_df[, c("term", "Estimate", "Q2.5", "Q97.5")]

ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
  geom_point() +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Fixed Effect", y = "Estimate",
       title = "Fixed Effects (Posterior Means & 95% CI)") +
  theme_minimal()



