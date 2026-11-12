
library(brms)
library(ggplot2)
library(patchwork)

mod <- readRDS("output/models/brm_nb_dt1_full_c_aint_ry.rds")
dt <- read.csv("output/dt1.csv")
dt$site <- factor(dt$site)
dt$year <- factor(dt$year)
dt$n_adu_log_c     <- as.numeric(scale(dt$n_adu_log,     scale = FALSE))
dt$mean_sett_log_c <- as.numeric(scale(dt$mean_sett_log, scale = FALSE))
dt$R_log10_c       <- as.numeric(scale(dt$R_log10,       scale = FALSE))
dt$Dc              <- as.numeric(scale(dt$D,             scale = FALSE))

preds <- c("n_adu_log_c", "Dc", "R_log10_c", "mean_sett_log_c")

labs_map <- c(
n_adu_log_c     = "Adult abundance (log, centered)",
Dc              = "Fractal dimension (centered)",
R_log10_c       = "Surface rugosity (log10, centered)",
mean_sett_log_c = "Average settlers (log, centered)"
)

# build aggregated prediction data for each pair
build_pair_data <- function(var_x, var_y, n_grid = 80) {
  grid_x <- seq(min(dt[[var_x]], na.rm = TRUE), max(dt[[var_x]], na.rm = TRUE), length.out = n_grid)
  grid_y <- seq(min(dt[[var_y]], na.rm = TRUE), max(dt[[var_y]], na.rm = TRUE), length.out = n_grid)
  pred_grid <- expand.grid(year  = levels(dt$year),
                           x_val = grid_x,
                           y_val = grid_y)
  for (p in preds)
    pred_grid[[p]] <- 0
  pred_grid[[var_x]] <- pred_grid$x_val
  pred_grid[[var_y]] <- pred_grid$y_val
  fit <- fitted(mod,
                newdata = pred_grid,
                re_formula = NA,
                summary = TRUE)
  pred_grid$estimate <- fit[, "Estimate"]
  pred_grid$ci_width <- fit[, "Q97.5"] - fit[, "Q2.5"]
  agg <- aggregate(
    cbind(estimate, ci_width) ~ get(var_x) + get(var_y),
    data = pred_grid,
    FUN = mean
  )
  names(agg)[1:2] <- c(var_x, var_y)
  agg
}
other_preds <- setdiff(preds, "n_adu_log_c")
pair_data <- lapply(other_preds, function(v) {
  list(var_y = v, data = build_pair_data("n_adu_log_c", v))
})
# compute shared color scale limits?
est_range <- range(unlist(lapply(pair_data, function(x)
  x$data$estimate)), na.rm = TRUE)
ci_range  <- range(unlist(lapply(pair_data, function(x)
  x$data$ci_width)), na.rm = TRUE)
# build panels using shared color scale
make_panel <- function(agg,
                       var_y,
                       fill_var,
                       fill_lab,
                       palette,
                       limits,
                       title) {
  ggplot(agg, aes_string(x = "n_adu_log_c", y = var_y, fill = fill_var)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(
      aes_string(z = fill_var),
      colour = "white",
      alpha = 0.4,
      linewidth = 0.3
    ) +
    geom_point(
      data = dt,
      aes_string(x = "n_adu_log_c", y = var_y),
      inherit.aes = FALSE,
      colour = "white",
      alpha = 0.25,
      size = 1
    ) +
    scale_fill_viridis_c(option = palette,
                         name = fill_lab,
                         limits = limits) +
    labs(x = labs_map[["n_adu_log_c"]], y = labs_map[[var_y]], title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text    = element_text(colour = "black"),
          panel.border = element_blank())
}
pair_rows <- lapply(pair_data, function(pd) {
  mean_panel <- make_panel(
    pd$data,
    pd$var_y,
    "estimate",
    "Juvenile\nabundance",
    "magma",
    est_range,
    "Posterior mean"
  )
  ci_panel   <- make_panel(pd$data,
                           pd$var_y,
                           "ci_width",
                           "95% CI width",
                           "viridis",
                           ci_range,
                           "Uncertainty")
  mean_panel + ci_panel
})
# stack rows and collect shared legends
i <- wrap_plots(pair_rows, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Interactions with adult abundance",
                  theme = theme(plot.title = element_text(face = "bold", size = 15))) &
  theme(legend.position = "right")
i
df_m <- effects_clean %>%
  filter(model == substr(mod_file, 15, 39)) %>%
  mutate(
    term_label = recode(
      term,
      "adults (log)" = "Adult abundance (log, centered)",
      "settlers (log)" = "Average settlers (log, centered)",
      "R log10" = "Surface rugosity (log10, centered)",
      "D" = "Fractal dimension (centered)",
      "year2023" = "Year (2023)",
      "adults (log):settlers (log)" = "Adult abundance × Average settlers",
      "adults (log):R log10" = "Adult abundance × Surface rugosity",
      "adults (log):D" = "Adult abundance × Fractal dimension",
    ),
    effect_type = if_else(grepl(":", term), "Interaction", "Main effect")
  )
mod_file <- "output/models/brm_nb_dt1_full_c_aint_ry.rds"
dt_file <- "output/dt1.csv"
df_m <- effects_clean %>%
  filter(model == substr(mod_file, 15, 39)) %>%
  mutate(
    term_label = recode(
      term,
      "adults (log)" = "Adult abundance (log,centered)",
      "settlers (log)" = "Average settlers (log,centered)",
      "R log10" = "Surface rugosity (log10,centered)",
      "D" = "Fractal dimension  (centered)",
      "year2023" = "Year (2023)",
      "adults (log):settlers (log)" = "Adult abundance × Average settlers",
      "adults (log):R log10" = "Adult abundance × Surface rugosity",
      "adults (log):D" = "Adult abundance × Fractal dimension",
    ),
    effect_type = if_else(grepl(":", term), "Interaction", "Main effect")
  )
df_m <- df_m %>%
  arrange(desc(effect_type == "Interaction"), desc(mean)) %>%
  mutate(term_label = factor(term_label, levels = rev(term_label))) %>%
  filter(type == "fixed")
e <- ggplot(df_m, aes(x = mean, y = term_label)) +
  geom_errorbarh(
    aes(xmin = l95, xmax = u95),
    height = 0,
    size = 0.8,
    alpha = 0.7
  ) +
  geom_errorbarh(aes(xmin = l80, xmax = u80), height = 0, size = 1.4) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  scale_x_continuous(limits = x_limits) +
  labs(x = "Model effect sizes", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.margin = margin(
      t = 12,
      r = 8,
      b = 8,
      l = 8
    )
  )
e <- e + xlim(-1, 1.5)
e


draws_long <- tibble(
  # Main effects
  `Adult abundance`                     = draws$b_n_adu_log_c,
  `Average settlers`                    = draws$b_mean_sett_log_c,
  `Surface rugosity`                    = draws$b_R_log10_c,
  `Fractal dimension`                   = draws$b_Dc,
  
  # Interactions with adult abundance
  `Adult abundance x Average settlers`  = draws$`b_n_adu_log_c:mean_sett_log_c`,
  `Adult abundance x Surface rugosity`  = draws$`b_n_adu_log_c:R_log10_c`,
  `Adult abundance x Fractal dimension` = draws$`b_n_adu_log_c:Dc`
) %>%
  pivot_longer(everything(), names_to = "predictor", values_to = "estimate") %>%
  mutate(
    type = if_else(str_detect(predictor, " x "), "Interaction", "Main effect"),
    predictor = factor(
      predictor,
      levels = c(
        "Adult abundance",
        "Average settlers",
        "Surface rugosity",
        "Fractal dimension",
        "Adult abundance x Average settlers",
        "Adult abundance x Surface rugosity",
        "Adult abundance x Fractal dimension"
      )
    )
  )


d <- ggplot(draws_long, aes(
  x = estimate,
  y = fct_rev(predictor),
  fill = type
)) +
  stat_halfeye(
    alpha = 0.7,
    point_interval = median_qi,
    .width = c(.5, .8, .95)
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(y = "Predictor", x = "Posterior estimate (log scale, negative binomial)", fill = NULL) +
  scale_fill_manual(values = c(
    "Main effect" = "#3A2767",
    "Interaction" = "#B19CD9"
  )) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")


# Need the centered adult-abundance values to set a sensible range
dt <- read.csv(dt_file)
dt$n_adu_log_c <- as.numeric(scale(dt$n_adu_log, scale = FALSE))

adu_grid <- seq(min(dt$n_adu_log_c, na.rm = TRUE),
                max(dt$n_adu_log_c, na.rm = TRUE),
                length.out = 60)

# For each interacting predictor X, the conditional slope at adult abundance a is:
# b_X + b_{n_adu_log_c:X} * a
build_slope <- function(main_col, int_col, label) {
  map_dfr(adu_grid, function(a) {
    slope <- draws[[main_col]] + draws[[int_col]] * a
    
    tibble(predictor   = label,
           n_adu_log_c = a,
           slope       = slope)
  })
}


slopes_long <- bind_rows(
  build_slope(
    "b_mean_sett_log_c",
    "b_n_adu_log_c:mean_sett_log_c",
    "Average settlers"
  ),
  build_slope("b_R_log10_c", "b_n_adu_log_c:R_log10_c", "Surface rugosity"),
  build_slope("b_Dc", "b_n_adu_log_c:Dc", "Fractal dimension")
)

# Summarise to median + 50/80/95% intervals
slopes_summary <- slopes_long %>%
  group_by(predictor, n_adu_log_c) %>%
  median_qi(slope, .width = c(.5, .8, .95)) %>%
  ungroup() %>%
  mutate(predictor = factor(
    predictor,
    levels = c("Average settlers", "Surface rugosity", "Fractal dimension")
  ))

cs <- ggplot(slopes_summary,
             aes(
               x = n_adu_log_c,
               y = slope,
               ymin = .lower,
               ymax = .upper
             )) +
  geom_hline(yintercept = 0,
             linetype = 2,
             colour = "grey40") +
  geom_lineribbon(aes(fill = factor(.width)), alpha = 0.5) +
  geom_line(linewidth = 0.7, colour = "#3A2767") +
  facet_wrap( ~ predictor, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "0.5"  = "#3A2767",
      "0.8"  = "#7A5BA8",
      "0.95" = "#B19CD9"
    ),
    name = "Credible interval"
  ) +
  labs(x = "Adult abundance (log, centered)", y = "Conditional slopes\nfor Juvenile abundance") +
  guides(fill = guide_legend(
    title.position = "left",
    title.hjust = 0.5,
    nrow = 1,
    byrow = TRUE
  )) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(
      t = 4,
      r = 0,
      b = 0,
      l = 0
    ),
    plot.margin = margin(
      t = 12,
      r = 8,
      b = 12,
      l = 8
    )
  )
cs

cs <- cs +
  geom_rug(
    data = dt,
    aes(x = n_adu_log_c),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.5
  )


build_slope <- function(main_col, int_col, label) {
  
  map_dfr(adu_grid, function(a) {
    
    tibble(
      draw        = seq_len(nrow(draws)),
      predictor   = label,
      n_adu_log_c = a,
      slope       = draws[[main_col]] + draws[[int_col]] * a
    )
  })
}

set.seed(123)

draw_subset <- slopes_long %>%
  distinct(predictor, draw) %>%
  group_by(predictor) %>%
  slice_sample(n = 1000)

cs +
  geom_line(
    data = semi_join(slopes_long,
                     draw_subset,
                     by = c("predictor", "draw")),
    aes(x = n_adu_log_c,
        y = slope,
        group = draw),
    inherit.aes = FALSE,
    colour = "black",
    alpha = 0.05,
    linewidth = 0.3
  )








library(tidyverse)
library(tidybayes)
library(ggdist)
library(patchwork)

plot_effect <- function(mod,
                        dt,
                        predictor,
                        xlab,
                        response = "n_juv") {
  
  dt$n_adu_log_c     <- as.numeric(scale(dt$n_adu_log,     scale = FALSE))
  dt$mean_sett_log_c <- as.numeric(scale(dt$mean_sett_log, scale = FALSE))
  dt$R_log10_c       <- as.numeric(scale(dt$R_log10,       scale = FALSE))
  dt$Dc              <- as.numeric(scale(dt$D,             scale = FALSE))
  
  ## Build prediction grid
  newdata <- tibble(
    n_adu_log_c      = mean(dt$n_adu_log_c),
    mean_sett_log_c  = mean(dt$mean_sett_log_c),
    R_log10_c        = mean(dt$R_log10_c),
    Dc               = mean(dt$Dc)
  )
  
  ## Replace focal predictor with a sequence
  newdata <- newdata[rep(1, 100), ]
  
  newdata[[predictor]] <- seq(
    min(dt[[predictor]], na.rm = TRUE),
    max(dt[[predictor]], na.rm = TRUE),
    length.out = 100
  )
  
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
      x = .data[[predictor]],
      y = .epred,
      ymin = .lower,
      ymax = .upper
    )
  ) +
    
    ## observed data
    geom_point(
      data = dt,
      aes(
        x = .data[[predictor]],
        y = .data[[response]]
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      alpha = .20,
      size = 1.2
    ) +
    
    ## credible intervals
    geom_lineribbon(
      aes(fill = factor(.width)),
      alpha = .50
    ) +
    
    ## median prediction
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
    ylim(0,110)+
    
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


p1 <- plot_effect(
  mod,
  dt,
  predictor = "n_adu_log_c",
  xlab = "Adult abundance (log, centred)"
)

p2 <- plot_effect(
  mod,
  dt,
  predictor = "R_log10_c",
  xlab = "Rugosity"
)

p3 <- plot_effect(
  mod,
  dt,
  predictor = "Dc",
  xlab = "Fractal dimension"
)

p2 +
  scale_color_grey()

12

p1
p2|p3





# Free the bottom panel from patchwork alignment with the wide y-axis labels above
p <- e / patchwork::free(cs) +
  plot_layout(heights = c(1, 1.35)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 18),
        plot.tag.position = c(0, 1))
p













#### maldives #######


############################################
#### code for Brambilla et. al analysis ####
############################################

## upload libraries
library(brms)
library(rstan)
library(ggplot2)
library(tidyverse)

## read data
df.tot<-read.csv("~/Library/CloudStorage/Dropbox/papers_conferences/Maldives recruitment study/Submission4_CoralReefs/versioning/recruitment_data_resubmit.csv", head = TRUE)

str(df.tot)
df.tot$srug <- scale(df.tot$sr) ## mean center surface rugosity 
df.tot$rack <- as.factor(df.tot$rack)   # factorize rack

#### fitting models ####

# Rstan settings

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### status, rugosity, depth ### 
set.seed(1234)
mod_stat_rug_depth <- brm(s.pres ~ alive.dead + srug + site + (1|rack),
                          data= df.tot, family = bernoulli(),
                          control = list(adapt_delta = 0.99),
                          iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_rug_depth)
plot(mod_stat_rug_depth) # chains check
#plot(marginal_effects(mod_stat_rug_depth)) # visualize effects

### status, depth ###
set.seed(1234)
mod_stat_depth <- brm(s.pres ~ alive.dead + site + (1|rack),
                      data= df.tot, family = bernoulli(),
                      control = list(adapt_delta = 0.99),
                      iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_depth)
plot(mod_stat_depth) # chains check
#plot(marginal_effects(mod_stat_depth)) # visualize effects

### rugosity, depth ###
set.seed(1234)
mod_rug_depth <- brm(s.pres ~ srug + site + (1|rack),
                     data= df.tot, family = bernoulli(),
                     control = list(adapt_delta = 0.99),
                     iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_depth)
plot(mod_rug_depth) # chains check
#plot(marginal_effects(mod_rug_depth)) # visualize effects

### rugosity, status ###
set.seed(1234)
mod_rug_stat <- brm(s.pres ~ srug + alive.dead + (1|rack),
                    data= df.tot, family = bernoulli(),
                    control = list(adapt_delta = 0.99),
                    iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_stat)
plot(mod_rug_stat) # chains check
#plot(marginal_effects(mod_rug_stat)) # visualize effect

### rugosity ###
set.seed(1234)
mod_rug <- brm(s.pres ~ srug + (1|rack),
               data= df.tot, family = bernoulli(),
               control = list(adapt_delta = 0.99),
               iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug)
plot(mod_rug) # chains check
#plot(marginal_effects(mod_rug)) # visualize effects

### status ###
set.seed(1234)
mod_stat <- brm(s.pres ~ alive.dead + (1|rack),
                data= df.tot, family = bernoulli(),
                control = list(adapt_delta = 0.99),
                iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat)
plot(mod_stat) # chains check
#plot(marginal_effects(mod_stat), ask = FALSE) # visualize effects

### depth ###
set.seed(1234)
mod_depth <- brm(s.pres ~ site + (1|rack),
                 data= df.tot, family = bernoulli(),
                 control = list(adapt_delta = 0.99),
                 iter = 20000, warmup = 1000, thin = 10)

summary(mod_depth)
plot(mod_depth) # chains check
#plot(marginal_effects(mod_depth), ask = FALSE) # visualize effects

### summary table ###
model.results <- data.frame(rbind(summary(mod_stat_rug_depth)$fixed,
                                  summary(mod_stat_depth)$fixed,
                                  summary(mod_rug_depth)$fixed,
                                  summary(mod_rug_stat)$fixed,
                                  summary(mod_rug)$fixed,
                                  summary(mod_stat)$fixed,
                                  summary(mod_depth)$fixed))

model.results <- data.frame(rbind(summary(mod_stat_rug_depth)$fixed,
                                  summary(mod_stat_depth)$fixed,
                                  summary(mod_rug_depth)$fixed,
                                  summary(mod_rug_stat)$fixed,
                                  summary(mod_rug)$fixed,
                                  summary(mod_stat)$fixed,
                                  summary(mod_depth)$fixed))

# compute and compare  WAIC and LOOic

WAIC(mod_stat_rug_depth, mod_stat_depth, mod_rug_depth, mod_rug_stat,mod_rug,
     mod_stat, mod_depth) # model 5 is the best fit

l1 <- loo(mod_stat_rug_depth)
l2 <- loo(mod_stat_depth)
l3 <- loo(mod_rug_depth)
l4 <- loo(mod_rug_stat)
l5 <- loo(mod_rug)
l6 <- loo(mod_stat)
l7 <- loo(mod_depth)

loo_compare(l1,l2,l3,l4,l5,l6,l7) # model 5 is the best fit

# create table
models <- c(rep("presence ~ alive.dead + srug + site + (1 | rackID)", 4),
            rep("presence ~ alive.dead + site + (1 | rackID)", 3),
            rep("presence ~ srug + site + (1 | rackID)", 3),
            rep("presence ~ srug + nubbin.status + (1 | rackID)", 3),
            rep("presence ~ srug + (1 | rackID)", 2),
            rep("presence ~ alive.dead + (1 | rackID)", 2),
            rep("presence ~ site + (1 | rackID)", 2))

model.results <- cbind(model = models, variable = rownames(model.results), model.results)

# write.csv(model.results,"results_summary.csv", row.names = FALSE, col.names = TRUE)

#### fig.3 ####

# add tile category based on combination of presence/absence of settlers 
# and status of coral nubbin
df.tot$rec.graph[df.tot$s.pres == 1 & df.tot$alive.dead == "alive"] <- "a"
df.tot$rec.graph[df.tot$s.pres == 1 & df.tot$alive.dead == "dead"] <- "b"
df.tot$rec.graph[df.tot$s.pres == 0 & df.tot$alive.dead == "alive"] <- "c"
df.tot$rec.graph[df.tot$s.pres == 0 & df.tot$alive.dead == "dead"] <-  "d"

# add sum of the tiles in each category and round sr at 3 digits
df.tot <- df.tot %>% group_by(rack) %>% group_by(alive.dead) %>% 
  mutate(spat.per = sum(s.pres), srr = round(sr, 3))

df.tot <- arrange(df.tot,sr,site)

# plot pie chart
ggplot(df.tot, aes(x=1, y= spat.per, fill=factor(rec.graph))) +
  geom_bar(stat="identity", position = position_fill()) +
  scale_fill_manual(values = viridis(9)[c(1,4,7,8)],
                    labels = c("  tiles with adult coral and settlers", 
                               "  tiles with dead coral and settlers", 
                               "  tiles with adult coral and no settlers", 
                               "  tiles with dead coral and no settlers"))+
  facet_wrap (site~srr, nrow = 2) +
  coord_polar(theta = "y", direction = -1)+
  ggtitle(" Tiles proportion in each rack\n")+
  theme_minimal()+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        title =  element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA,), # bg of the plot
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel


#### fig.4 ####

## fig 4a - marginal effect
m <- plot(marginal_effects(mod_rug, col = c("black","red")))
fig4a <- m[[1]] +
  geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
  theme_minimal()+ ylim(c(0,1))+
  scale_color_grey() +
  scale_fill_grey() +
  labs(x = "\n standardized rugosity", y = "presence probability \n") +
  #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
  #           position = position_jitter(width = 0.01)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.position ="top", 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel

# fig 4b - posterior distribution of the slopes of mod_rug
posterior_slope <- as.array(mod_rug)
dim(posterior_slope)
dimnames(posterior_slope) 
p <-plot(posterior_slope, pars = "b_srug")

library(bayesplot)
d <- mcmc_areas(
  posterior_slope, 
  pars = "b_srug",
  prob = .95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) 

d + 
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_discrete(breaks="b_srug",
                   labels="surface \n rugosity")+
  labs(x = "\n slope", y = " frequency\n")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        legend.text = element_text(size = 12),legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.position ="top", 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel

# fig 4c - marginal effect of status of the nubbin in mod_stat

s <- plot(marginal_effects(mod_stat))

fig4c <- s[[1]] +
  #geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
  #theme_minimal()+ ylim(c(0,1))+
  #scale_color_grey() +
  #scale_fill_grey() +
  labs(x = "\n adult coral status", y = "presence probability \n") +
  #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
  #           position = position_jitter(width = 0.01)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.position =c(0.8,0.6), 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_blank())) # get rid of legend panel)) # get rid of legend panel


### ICRS presentation

ce <- conditional_effects(mod_rug)
plot_df <- ce$srug

ggplot(plot_df, aes(x = srug, y = estimate__)) +
geom_ribbon(
aes(ymin = lower__, ymax = upper__),
fill = "#B19CD9",
alpha = .8
) +
geom_line(
colour = "#3A2767",
linewidth = 2
) +
  geom_hline(yintercept = c(0, 1), colour = "grey90", size = .4)+
geom_count(
data = df.tot,
aes(x = srug, y = s.pres),
inherit.aes = FALSE,
fill = "#3A2767",
col = "#B19CD9",
alpha = 0.5
) +
coord_cartesian(ylim = c(0, 1)) +
labs(
x = "Standardized rugosity",
y = "Settler presence probability",
size = "# tiles"
) +
theme_minimal(base_size = 16)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




df.tot$alive.dead <- as.factor(df.tot$alive.dead)

ce_stat <- conditional_effects(mod_stat)

newdat <- data.frame(
  alive.dead = levels(df.tot$alive.dead)
)
pred <- fitted(
  mod_stat,
  newdata = newdat,
  summary = TRUE,
  scale = "response"
)
newdat <- data.frame(
  alive.dead = levels(df.tot$alive.dead)
)
pred <- fitted(
  mod_stat,
  newdata = newdat,
  re_formula = NA,      # <- ignore random effects
  scale = "response"
)
plot_df <- cbind(newdat, as.data.frame(pred))

ggplot(plot_df,
aes(x = alive.dead, y = Estimate)) +
geom_errorbar(
aes(ymin = Q2.5, ymax = Q97.5),
colour = "#3A2767",
width = 0.15,
linewidth = 0.8
) +
geom_point(
shape = 21,
size = 4,
fill = "#B19CD9",
colour = "#3A2767",
stroke = 1
) +
labs(
x = "Status",
y = "Settler presence probability"
) +
theme_classic(base_size = 14)
ggplot(plot_df,
aes(x = alive.dead, y = Estimate)) +
geom_errorbar(
aes(ymin = Q2.5, ymax = Q97.5),
colour = "#3A2767",
width = 0.15,
linewidth = 0.8
) +
geom_point(
shape = 21,
size = 4,
fill = "#B19CD9",
colour = "#3A2767",
stroke = 1
) +
labs(
x = "Adult status",
y = "Settler presence probability"
) +
theme_classic(base_size = 16)

