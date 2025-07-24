### merge counts and complexity 

library(dplyr)
library(ggplot2)
library(sjPlot)
library(mgcv)
library(lme4)
library(glmmTMB)
library(DHARMa)
# for soap smoother
library(sf)
library(sp)
library(rgeos)
library(terra)
library(broom.mixed)

###
dt1 <- read.csv("dt1.csv")
dt1$site <- as.factor(dt1$site)
dt1$year <- as.factor(dt1$year)

dt1 %>% arrange(desc(H))
dt1 <- dt1[dt1$H != max(dt1$H),]

dt1 %>% group_by(site,year)%>%mutate(av_R = mean(R)) %>% select(site, year,av_R) %>% unique() %>% arrange(year,av_R)
dt1 %>% group_by(site,year)%>%mutate(av_adu = mean(n_adu)) %>% select(site, year,av_adu) %>% unique() %>% arrange(year,av_adu)
dt1 %>% group_by(site,year)%>%mutate(av_juv = mean(n_juv)) %>% select(site, year,av_juv) %>% unique() %>% arrange(year,av_juv)

summary_data <- dt1 %>%
  group_by(site) %>%
  summarise(mean_count = mean(n_juv),
            var_count = var(n_juv))

plot(summary_data$mean_count, summary_data$var_count,
     xlab = "Mean", ylab = "Variance", main = "Mean-Variance Plot across sites")



#fit poisson
poisson_model <- glmer(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + year +
                         (1 | site), 
                       data = dt1,
                       family = poisson)

resid_pearson <- residuals(poisson_model, type = "pearson")
df_resid <- nrow(dt1) - length(fixef(poisson_model))
dispersion_ratio <- sum(resid_pearson^2) / df_resid
dispersion_ratio #overdispersed

# fit negative binomial
nb_model1 <- glmmTMB(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + year +
                      (1 | site), 
                    data = dt1,
                    family = nbinom1)
nb_model2 <- glmmTMB(n_juv ~ n_adu_log + mean_sett_s + R_log10 + D_s + year +
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

gamm1 <- gam(data = dt1, n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year, correlation = corGaus(form = ~ ys + xs),
             method="REML", family = ziP())
summary(gamm1)
tab_model(gamm1)
gam.check(gamm1, old.style = TRUE)
plot(residuals(gamm1))

gamm1_noyear <- gam(data = dt1, n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s, correlation = corGaus(form = ~ ys + xs),
             method="REML", family = ziP())
summary(gamm1_noyear)
tab_model(gamm1_noyear)
gam.check(gamm1_noyear, old.style = TRUE)
plot(residuals(gamm1_noyear))

### build soap smoother ####

# read island shp and creat soap area
ann<- st_read("data_raw/ann_shifted_manual.shp")
island_sf <- st_read("data_raw/island/dyiigurra_boundary.shp")
reef_sites <- st_read("data_raw/island/site_centroids.shp")
island_sf <- st_transform(island_sf, crs(ann))
reef_sites  <- st_transform(reef_sites , crs(ann))
# plot(st_geometry(island_sf), col = 'grey', main = "Dyiigurra")
# points(reef_sites)

reef_buffer <- st_buffer(island_sf, dist = 1500)  # 2 km buffer
reef_union <- st_union(reef_buffer)
island_union <- st_union(island_sf)

reef_buffer_grid <- st_buffer(island_sf, dist = 1400)  # 1.4 km buffer
reef_buffer_grid_is <- st_buffer(island_sf, dist = 50)  # 10m buffer
reef_buffer_grid_union <- st_union(reef_buffer_grid_is)
reef_union_grid <- st_union(reef_buffer_grid)

soap_domain <- st_difference(reef_union, island_union)

grid_domain <- st_difference(reef_union_grid, reef_buffer_grid_union)


get_xy <- function(polygon) {
  # Take the outer ring only
  coords <- st_coordinates(polygon)[, 1:2]
  list(x = coords[,1], y = coords[,2])
}

# Convert outer boundary (soap_domain)
soap_boundary <- get_xy(soap_domain)
outer_bnd <- get_xy(reef_union)
# Convert each island polygon to hole
island_holes <- lapply(1:nrow(island_sf), function(i) {
  poly <- island_sf[i, ]
  get_xy(poly)
})

# Final bnd list: outer first, then islands
bnd <- c(list(outer_bnd), island_holes)


grid <- expand.grid(
  x = seq(min(dt1$x0)-1500, max(dt1$x0)+1500, length = 50),
  y = seq(min(dt1$y0)-1900, max(dt1$y0)+1500, length = 50)
)

grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = crs(ann))
grid_in_domain <- st_intersection(grid_sf, grid_domain)

plot(soap_domain, col = "grey50")
points(grid_in_domain, col = "blue")
knots <- st_coordinates(grid_in_domain)
colnames(knots) <- c("x","y")


####### 1x1m scale ##########


# soap c0 ####

dt1$x <- dt1$x0
dt1$y <- dt1$y0


x <- dt1$x
y <- dt1$y

valid_idx <- is.finite(x) & is.finite(y)
dt1 <- dt1[valid_idx, ]

inside_data <- inSide(bnd, x, y)
table(inside_data)
dt1$inside <- inside_data
dt1 <- dt1[inside_data, ]


soap_c0 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = nb())


summary(soap_c0)
tab_model(soap_c0)

soap_c0_ziP2<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                   ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                   data = dt1, method = "REML",
                 knots = list(x = knots[,1], y = knots[,2]),
                 family = ziP())

soap_c0_ziP2_int<- mgcv::gam(n_juv ~ year*n_adu_log + year*R_log10 + year*D_s+ year*mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                         data = dt1, method = "REML",
                         knots = list(x = knots[,1], y = knots[,2]),
                         family = ziP())

soap_c0_ziP2_19<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                         ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                         data = dt1[dt1$year == "2019",], method = "REML",
                         knots = list(x = knots[,1], y = knots[,2]),
                         family = ziP())

soap_c0_ziP2_22<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                         ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                         data = dt1[dt1$year =="2022",], method = "REML",
                         knots = list(x = knots[,1], y = knots[,2]),
                         family = ziP())


soap_c0_ziP <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                   data = dt1, method = "REML",
                   knots = list(x = knots[,1], y = knots[,2]),
                   family = ziP())

soap_c0_22 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2022",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())

soap_c0_19 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2019",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())

summary(soap_c0_19)
summary(soap_c0_22)
tab_model(soap_c0, soap_c0_19, soap_c0_22)


saveRDS(soap_c0,"outputs/soap_c0.rds")
saveRDS(soap_c0_ziP,"outputs/soap_c0_ziP.rds")
saveRDS(soap_c0_19,"outputs/soap_c0_19.rds")
saveRDS(soap_c0_22,"outputs/soap_c0_22.rds")

saveRDS(soap_c0,"outputs/soap_c0.rds")
saveRDS(soap_c0_ziP2,"outputs/soap_c0_ziP2.rds")
saveRDS(soap_c0_ziP2_19,"outputs/soap_c0_ziP2_19.rds")
saveRDS(soap_c0_ziP2_22,"outputs/soap_c0_ziP2_22.rds")

vis.gam(soap_c0, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_19, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_22, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_ziP2, view = c("x", "y"), plot.type = "contour", asp = 1)


dt1$x <- dt1$xs
dt1$y <- dt1$ys

soap_cs <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = nb())
summary(soap_cs)
tab_model(soap_cs)

soap_cs_ziP <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                 data = dt1, method = "REML",
                 knots = list(x = knots[,1], y = knots[,2]),
                 family = ziP())


soap_cs_22 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2022",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())


soap_cs_19 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2019",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())






saveRDS(soap_cs,"outputs/soap_cs.rds")
saveRDS(soap_cs_ziP,"outputs/soap_cs_ziP.rds")
saveRDS(soap_cs_19,"outputs/soap_cs_19.rds")
saveRDS(soap_cs_22,"outputs/soap_cs_22.rds")
ggplot() +
  geom_sf(data = st_as_sf(soap_domain), fill = NA, color = "black") +
  geom_point(data = dt1, aes(x = x, y = y, color = inside)) +
  scale_color_manual(values = c("red", "green")) +
  theme_minimal()


summary(soap1_k10)
tab_model(soap1)
gam.check(soap1, old.style = TRUE)
plot(residuals(soap1))
gratia::draw(soap1)

### plot predictions

grid <- expand.grid(
  x = seq(min(dt1$x) - 2000, max(dt1$x) + 2000, length = 100),
  y = seq(min(dt1$y) - 2000, max(dt1$y) + 2000, length = 100)
)

grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(ann))
grid_in_domain <- st_intersection(grid_sf, soap_domain)
grid_coords <- st_coordinates(grid_in_domain)

# Create prediction data frame
pred_data <- as.data.frame(grid_coords)
colnames(pred_data) <- c("x", "y")

# Add mean or representative values for all other model predictors
pred_data$n_adu_log    <- mean(dt1$n_adu_log, na.rm = TRUE)
pred_data$R_log10      <- mean(dt1$R_log10, na.rm = TRUE)
pred_data$D_s          <- mean(dt1$D_s, na.rm = TRUE)
pred_data$mean_sett_s  <- mean(dt1$mean_sett_s, na.rm = TRUE)
pred_data$year         <- as.factor("2019")

pred_data$fit <- predict(soap1, newdata = pred_data, type = "response")

ggplot(pred_data, aes(x,y, fill = fit)) +
  geom_raster(interpolate = TRUE)+
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  labs(title = "Predicted juvenile counts", fill = "Prediction")


ggplot(pred_data, aes(x = x, y = y, z = fit)) +
  geom_contour_filled() +
  coord_equal() +
  theme_minimal()

### soap 2 ####
dt1$x <- dt1$x0
dt1$y <- dt1$y0

grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = crs(ann))
grid_in_domain <- st_intersection(grid_sf, grid_domain)

## create knot grid
grid <- expand.grid(
  x = seq(min(dt1$x0)-2000, max(dt1$x0)+2000, length = 70),
  y = seq(min(dt1$y0)-2000, max(dt1$y0)+2000, length = 70)
)


plot(soap_domain)
points(grid_in_domain, col = "blue")
knots <- st_coordinates(grid_in_domain)
colnames(knots) <- c("x","y")

x <- dt1$x
y <- dt1$y
inside_data <- inSide(bnd, x, y)
table(inside_data)
dt1$inside <- inside_data
dt1 <- dt1[inside_data, ]





soap2 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = ziP())

saveRDS(soap2, file = "outputs/soap2.rds")

soap2_19 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == 2019,], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = ziP())

saveRDS(soap2_19, file = "outputs/soap2_19.rds")

soap2_22 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == 2022,], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = ziP())

saveRDS(soap2_22, file = "outputs/soap2_22.rds")

summary(soap2)
summary(soap2_19)
summary(soap2_22)

ggplot() +
  geom_sf(data = st_as_sf(soap_domain), fill = NA, color = "black") +
  geom_point(data = dt1, aes(x = x, y = y, color = inside)) +
  scale_color_manual(values = c("red", "green")) +
  theme_minimal()

saveRDS(soap2_22, file = "outputs/soap2_22.rds")

summary(soap2)
tab_model(soap2)
gam.check(soap2, old.style = TRUE)
plot(residuals(soap2))
gratia::draw(soap2)

### plot predictions

grid <- expand.grid(
  x = seq(min(dt1$x) - 2000, max(dt1$x) + 2000, length = 80),
  y = seq(min(dt1$y) - 2000, max(dt1$y) + 2000, length = 80)
)

grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(ann))
grid_in_domain <- st_intersection(grid_sf, soap_domain)
grid_coords <- st_coordinates(grid_in_domain)

# Create prediction data frame
pred_data <- as.data.frame(grid_coords)
colnames(pred_data) <- c("x", "y")

# Add mean or representative values for all other model predictors
pred_data$n_adu_log    <- mean(dt1$n_adu_log, na.rm = TRUE)
pred_data$R_log10      <- mean(dt1$R_log10, na.rm = TRUE)
pred_data$D_s          <- mean(dt1$D_s, na.rm = TRUE)
pred_data$mean_sett_s  <- mean(dt1$mean_sett_s, na.rm = TRUE)
pred_data$year         <- as.factor("2019")

pred_data$fit <- predict(soap2, newdata = pred_data, type = "response")

ggplot(pred_data, aes(x,y, fill = fit)) +
  geom_raster(interpolate = TRUE)+
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  labs(title = "Predicted juvenile counts", fill = "Prediction")


ggplot(pred_data, aes(x = x, y = y, z = fit)) +
  geom_contour_filled() +
  coord_equal() +
  theme_minimal()


#### soap 3 ###


dt1$x <- dt1$x0
dt1$y <- dt1$y0


soap3 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(site, bs = "re") +
               s(x, y, bs = "so", xt = list(bnd = bnd)),
data = dt1, method = "REML",
knots = list(x = knots[,1], y = knots[,2]),
family = ziP())

saveRDS(soap3, file = "outputs/soap3.rds")

soap4 <- gam(n_juv ~ n_adu_log*R_log10 + n_adu_log*D_s+ mean_sett_s + year+ s(site, bs = "re"),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = ziplss())

summary(soap4)

soap5 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = nb())
summary(soap5)
tab_model(soap5)

soap5_ziP <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = ziP())


soap5_22 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                 data = dt1[dt1$year == "2022",], method = "REML",
                 knots = list(x = knots[,1], y = knots[,2]),
                 family = nb())


soap5_19 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                 data = dt1[dt1$year == "2019",], method = "REML",
                 knots = list(x = knots[,1], y = knots[,2]),
                 family = nb())


summary(soap5_19)
summary(soap5_22)
tab_model(soap5, soap5_19, soap5_22)



dt1$x <- dt1$xs
dt1$y <- dt1$ys

soap6 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = nb())
summary(soap5)
tab_model(soap5)

soap6_ziP <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                 data = dt1, method = "REML",
                 knots = list(x = knots[,1], y = knots[,2]),
                 family = ziP())


soap6_22 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2022",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())


soap6_19 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                data = dt1[dt1$year == "2019",], method = "REML",
                knots = list(x = knots[,1], y = knots[,2]),
                family = nb())





saveRDS(soap5,"outputs/soap5.rds")
saveRDS(soap5_ziP,"outputs/soap5_ziP.rds")
saveRDS(soap5_19,"outputs/soap5_19.rds")
saveRDS(soap5_22,"outputs/soap5_22.rds")

saveRDS(soap6,"outputs/soap6.rds")
saveRDS(soap6_ziP,"outputs/soap6_ziP.rds")
saveRDS(soap6_19,"outputs/soap6_19.rds")
saveRDS(soap6_22,"outputs/soap6_22.rds")

dt1$x <- dt1$xs
dt1$y <- dt1$ys

soap6 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
             data = dt1, method = "REML",
             knots = list(x = knots[,1], y = knots[,2]),
             family = nb())
summary(soap6)
tab_model(soap6)
vis.gam(soap_c0, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_19, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_22, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(soap_c0_ziP, view = c("x", "y"), plot.type = "contour", asp = 1)


plot(soap5, scheme = 2, asp = 1, pages = 1)
plot(soap5_19, scheme = 2, asp = 1, pages = 1)
plot(soap5_22, scheme = 2, asp = 1, pages = 1)


saveRDS(soap6,"outputs/soap6.rds")



### plot predictions

grid <- expand.grid(
  x = seq(min(dt1$x) - 2000, max(dt1$x) + 2000, length = 80),
  y = seq(min(dt1$y) - 2000, max(dt1$y) + 2000, length = 80)
)

grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(ann))
grid_in_domain <- st_intersection(grid_sf, soap_domain)
grid_coords <- st_coordinates(grid_in_domain)

# Create prediction data frame
pred_data <- as.data.frame(grid_coords)
colnames(pred_data) <- c("x", "y")

# Add mean or representative values for all other model predictors
pred_data$n_adu_log    <- mean(dt1$n_adu_log, na.rm = TRUE)
pred_data$R_log10      <- mean(dt1$R_log10, na.rm = TRUE)
pred_data$D_s          <- mean(dt1$D_s, na.rm = TRUE)
pred_data$mean_sett_s  <- mean(dt1$mean_sett_s, na.rm = TRUE)
pred_data$year         <- as.factor("2019")

pred_data$fit5 <- predict(soap5, newdata = pred_data, type = "response")

ggplot(pred_data, aes(x,y, fill = fit5)) +
  geom_raster(interpolate = TRUE)+
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  labs(title = "Predicted juvenile counts soap 5", fill = "Prediction")

pred_data$fit6 <- predict(soap6, newdata = pred_data, type = "response")

ggplot(pred_data, aes(x,y, fill = fit6)) +
  geom_raster(interpolate = TRUE)+
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  labs(title = "Predicted juvenile counts soap 6", fill = "Prediction")


# Plot histogram
hist(dt1$n_juv, breaks = 50, main = "Response Variable", xlab = "n_juv")

# Calculate proportion of zeros
mean(dt1$n_juv == 0)


# Fit standard models
m_pois <- glm(n_juv ~ n_adu_log +R_log10 + D_s+ mean_sett_s + year + site, data = dt1, family = poisson)
m_nb <- MASS::glm.nb(n_juv ~ n_adu_log +R_log10 + D_s+ mean_sett_s + year + site, data = dt1)

# Use DHARMa to test for zero inflation
library(DHARMa)
sim_pois <- simulateResiduals(m_pois)
testZeroInflation(sim_pois)

sim_nb <- simulateResiduals(m_nb)
testZeroInflation(sim_nb)


c0_22 <- readRDS("outputs/soap_c0_ziP2_22.rds")
c0_19 <- readRDS("outputs/soap_c0_ziP2_19.rds")
c0 <- readRDS("outputs/soap_c0_ziP2.rds")

c0
c0_19
c0_22

tab_model(c0_19, c0_22)
vis.gam(c0_19, view = c("x", "y"), plot.type = "contour", asp = 1)
vis.gam(c0_22, view = c("x", "y"), plot.type = "contour", asp = 1)
gam.check(c0_19, old.style = TRUE)
plot(c0_19, asp = 1)
points(dt1[dt1$year == 2019,]$x0,dt1[dt1$year == 2019,]$y0)
plot(c0_22, asp = 1)
points(dt1[dt1$year == 2022,]$x0,dt1[dt1$year == 2022,]$y0)

cs <- readRDS("outputs/soap6.rds")


c0_nb22 <- readRDS("outputs/soap5_22.rds")
c0_nb19 <- readRDS("outputs/soap5_19.rds")
c0_nb <- readRDS("outputs/soap5.rds")

plot(c0_nb19, asp = 1)
points(dt1[dt1$year == 2019,]$x0,dt1[dt1$year == 2019,]$y0)
plot(c0_nb22, asp = 1)
points(dt1[dt1$year == 2022,]$x0,dt1[dt1$year == 2022,]$y0)

tab_model(c0_nb19, c0_nb22)




####### 2x2m scale ##########

dt2 <- read.csv("dt2.csv")
dt2$site <- as.factor(dt2$site)
dt2$year <- as.factor(dt2$year)

dt2 %>% arrange(desc(H))
dt2 <- dt2[dt2$H != max(dt2$H),]

dt2 %>% group_by(site,year)%>%mutate(av_R = mean(R)) %>% select(site, year,av_R) %>% unique() %>% arrange(year,av_R)
dt2 %>% group_by(site,year)%>%mutate(av_adu = mean(n_adu)) %>% select(site, year,av_adu) %>% unique() %>% arrange(year,av_adu)
dt2 %>% group_by(site,year)%>%mutate(av_juv = mean(n_juv)) %>% select(site, year,av_juv) %>% unique() %>% arrange(year,av_juv)

summary_data <- dt2 %>%
  group_by(site) %>%
  summarise(mean_count = mean(n_juv),
            var_count = var(n_juv))

plot(summary_data$mean_count, summary_data$var_count,
     xlab = "Mean", ylab = "Variance", main = "Mean-Variance Plot across sites")

plot(dt2$n_juv)


dt2 <- dt2 %>% mutate(n_adu_log = log(n_adu),
                      R_log10 = log10(R),
                      D_s = scale(D),
                      mean_sett_s = scale(mean_sett))

# soap c0 ####

dt2$x <- dt2$x0
dt2$y <- dt2$y0


x <- dt2$x
y <- dt2$y

valid_idx <- is.finite(x) & is.finite(y)
dt2 <- dt2[valid_idx, ]

inside_data <- inSide(bnd, x, y)
table(inside_data)
dt2$inside <- inside_data
dt2 <- dt2[inside_data, ]


soap_c0 <- gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
               data = dt2, method = "REML",
               knots = list(x = knots[,1], y = knots[,2]),
               family = nb())


summary(soap_c0)
tab_model(soap_c0)

soap_c0_ziP2<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                         ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + year+ s(x, y, bs = "so", xt = list(bnd = bnd)),
                         data = dt2, method = "REML",
                         knots = list(x = knots[,1], y = knots[,2]),
                         family = ziP())

soap_c0_ziP2_int<- mgcv::gam(n_juv ~ year*n_adu_log + year*R_log10 + year*D_s+ year*mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                             data = dt2, method = "REML",
                             knots = list(x = knots[,1], y = knots[,2]),
                             family = ziP())

soap_c0_ziP2_19<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                            ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                            data = dt2[dt2$year == "2019",], method = "REML",
                            knots = list(x = knots[,1], y = knots[,2]),
                            family = ziP())

soap_c0_ziP2_22<- mgcv::gam(n_juv ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                            ziplp = ~ n_adu_log + R_log10 + D_s+ mean_sett_s + s(x, y, bs = "so", xt = list(bnd = bnd)),
                            data = dt2[dt2$year =="2022",], method = "REML",
                            knots = list(x = knots[,1], y = knots[,2]),
                            family = ziP())


     
