### merge counts and complexity 

library(dplyr)
library(ggplot2)
library(sjPlot)
library(mgcv)
library(lme4)
library(glmmTMB)
library(DHARMa)


###

dt1$site <- as.factor(dt1$site)

summary_data <- dt1 %>%
  group_by(site) %>%
  summarise(mean_count = mean(n_juv),
            var_count = var(n_juv))

plot(summary_data$mean_count, summary_data$var_count,
     xlab = "Mean", ylab = "Variance", main = "Mean-Variance Plot across sites")

#fit poisson
poisson_model <- glmer(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + 
                         (1 | site), 
                       data = dt1,
                       family = poisson)

resid_pearson <- residuals(poisson_model, type = "pearson")
df_resid <- nrow(dt1) - length(fixef(poisson_model))
dispersion_ratio <- sum(resid_pearson^2) / df_resid
dispersion_ratio #overdispersed

# fit negative binomial
nb_model1 <- glmmTMB(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + 
                      (1 | site), 
                    data = dt1,
                    family = nbinom1)
nb_model2 <- glmmTMB(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + 
                      (1 | site), 
                    data = dt1,
                    family = nbinom2)

res_nb1 <- simulateResiduals(nb_model1)
plot(res_nb1)
testDispersion(res_nb1)

res_nb2 <- simulateResiduals(nb_model2)
plot(res_nb2)
testDispersion(res_nb2)

AIC(nb_model1,nb_model2) # model 2 has better plots and lower AIC
summary(nb_model2)
tab_model(nb_model2)

ranef(nb_model2)
ran <- ranef(nb_model2)$cond$site
ran$site <- rownames(ran)

plot_model(nb_model2, type = "re")
ggplot(ran, aes(x = site, y = `(Intercept)`)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(y = "Random Intercept Effect (Site)", x = "Site")

# try with different transforms
dt1$n_adu_c <- dt1$n_adu - mean(dt1$n_adu)
dt1$mean_sett_c <- dt1$mean_sett - mean(dt1$mean_sett)
dt1$R_c <- dt1$R - mean(dt1$R)
dt1$D_c <- dt1$D - mean(dt1$D)

nb_model2_c <- glmmTMB(n_juv ~ n_adu_c + mean_sett_c + R_c + D_c + 
                       (1 | site), 
                     data = dt1,
                     family = nbinom2)

res_nb2_c <- simulateResiduals(nb_model2_c)
plot(res_nb2_c)
testDispersion(res_nb2_c)

summary(nb_model2_c)
tab_model(nb_model2_c)

ranef(nb_model2_c)
ran <- ranef(nb_model2_c)$cond$site
ran$site <- rownames(ran)
plot_model(nb_model2_c, type = "re")
ggplot(ran, aes(x = site, y = `(Intercept)`)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(y = "Random Intercept Effect (Site)", x = "Site")

# consistent

### start gamming 

gamm1 <- gam(data = dt1, n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s, correlation = corGaus(form = ~ ys + xs),
             method="REML", family = ziP())
summary(gamm1)
tab_model(gamm1)
gam.check(gamm1, old.style = TRUE)
plot(residuals(gamm1))
