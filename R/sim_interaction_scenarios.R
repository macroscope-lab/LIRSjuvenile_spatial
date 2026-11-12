# Simulate z from x and y under three scenarios, then plot all three heatmaps in one row:
#   1. Positive interaction  (b3 > 0)
#   2. Negative interaction  (b3 < 0)
#   3. Additive effects      (b3 = 0, no interaction)
#
# z = b0 + b1*x + b2*y + b3*x*y + noise
# Needs only ggplot2.

library(ggplot2)

set.seed(42)

# 1. Settings ---------------------------------------------------------------
n  <- 500     # number of simulated points per scenario
b0 <- 0       # intercept
b1 <- 1       # effect of x
b2 <- 1       # effect of y
noise_sd <- 0.3

scenarios <- c(
  "Positive interaction" =  1.5,
  "Negative interaction" = -1.5,
  "Additive (no interaction)" = 0
)

# 2. Simulate, fit and predict for each scenario ----------------------------
# x and y are centred on 0 so that the interaction is easy to see.
grid <- expand.grid(
  x = seq(-1, 1, length.out = 100),
  y = seq(-1, 1, length.out = 100)
)

sim_one <- function(name, b3) {
  x <- runif(n, -1, 1)
  y <- runif(n, -1, 1)
  z <- b0 + b1 * x + b2 * y + b3 * x * y + rnorm(n, sd = noise_sd)
  dat <- data.frame(x, y, z)

  # Fit the model WITH an interaction term in every case.
  # In the additive scenario the interaction estimate should come out near zero.
  fit <- lm(z ~ x * y, data = dat)
  cat("\n====", name, "====\n")
  print(summary(fit)$coefficients)

  g <- grid
  g$z_pred   <- predict(fit, newdata = g)
  g$scenario <- name
  dat$scenario <- name

  list(points = dat, grid = g)
}

results <- Map(sim_one, names(scenarios), scenarios)

points_all <- do.call(rbind, lapply(results, `[[`, "points"))
grid_all   <- do.call(rbind, lapply(results, `[[`, "grid"))

# Keep the panels in the order given above
lvls <- names(scenarios)
points_all$scenario <- factor(points_all$scenario, levels = lvls)
grid_all$scenario   <- factor(grid_all$scenario,   levels = lvls)

# 3. Three heatmaps in one row, sharing one colour scale ---------------------
# A shared scale lets you compare the scenarios directly.
p <- ggplot() +
  geom_raster(data = grid_all, aes(x = x, y = y, fill = z_pred)) +
  geom_point(data = points_all, aes(x = x, y = y), size = 0.5, alpha = 0.35) +
  scale_fill_viridis_c(name = "Predicted z") +
  facet_wrap(~ scenario, nrow = 1) +
  coord_fixed() +
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

print(p)
ggsave("heatmaps_interaction_scenarios.png", p, width = 13, height = 4.5, dpi = 300)
