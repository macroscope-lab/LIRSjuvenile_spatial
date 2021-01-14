#### preliminary analysis ####
library(mgcv)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(rgdal)
library(sp)
library(sf)
library(sjPlot)
library(sgeostat)

## read tables 

tab05 <- read.csv("data_tables/tab05_all.csv", head = TRUE)
tab1 <- read.csv("data_tables/tab1_all.csv", head = TRUE)
tab2 <- read.csv("data_tables/tab2_all.csv", head = TRUE)

summary(tab05)
dim(tab05)

tab2$rec<- as.factor(tab2$rec)
tab1$rec<- as.factor(tab1$rec)
tab2$rec<- as.factor(tab2$rec)

# check D and R values
plot(data = tab05, D~R)
tab05[tab05$R == max(tab05$R),]
tab05 <- tab05[-which(tab05$R == max(tab05$R)),]
plot(data = tab05, D~R) #better

#remove D<2?

plot(data = tab1, D~R)
tab1[tab1$R == max(tab1$R),]
tab1 <- tab1[-which(tab1$R == max(tab1$R)),]
plot(data = tab1, D~R)

plot(data = tab2, D~R)

## N3 map example
N3 <- raster("data_raw/records/N319.tif")
N3 <- projectRaster(N3,
                    crs = "+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ann <- readOGR("data_raw/annotations/Export_Output_all.shp")
ann_old <- subset(ann, ann$Status == "alive")
ann_new <- subset(ann, ann$Status == "new")

# DEM palette
ncols <- 100
newcol <- colorRampPalette(c("#8C510A", "#BF812D",
                             "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"))
cols <- newcol(ncols)  #apply the function to get 100 colours

Ccoord2 <- tab2[,c("xC", "yC")]
c2 <- SpatialPoints(Ccoord2)
Ccoord1 <- tab1[,c("xC", "yC")]
c1 <- SpatialPoints(Ccoord1)
Ccoord2 <- tab2[,c("xC", "yC")]
c2 <- SpatialPoints(Ccoord2)

#### grid parameters
g2 <- list()
for(i in seq_along(c2)) {
  g2[i] <- gBuffer(c2[i], width=1, quadsegs=1, capStyle="SQUARE")
}
G2 <- do.call(bind, g2)
crs(G2) <- crs(ann)

g1 <- list()
for(i in seq_along(c1)) {
  g1[i] <- gBuffer(c1[i], width=0.5, quadsegs=1, capStyle="SQUARE")
}
G1 <- do.call(bind, g1)
crs(G1) <- crs(ann)

g2 <- list()
for(i in seq_along(c2)) {
  g2[i] <- gBuffer(c2[i], width=0.25, quadsegs=1, capStyle="SQUARE")
}
G2 <- do.call(bind, g2)
crs(G2) <- crs(ann)

### get N3 center ###
n3 <- data.frame(lon = 615.239532, lat = 5865.968985)
center <- st_as_sf(x = n3, 
                   coords = c("lon", "lat"),
                   crs = crs(N3))
center_sp<- as(center, "Spatial")
N3c <- gBuffer(center_sp, width=7.3)
circle_raster <- mask(N3, N3c)

png("figs/CountMap2.png", width = 16, height = 14, unit = 'cm', res = 300)
plot(circle_raster, asp=1, col=cols, axes=FALSE, box=FALSE)
plot(ann_adults, add = TRUE, col = "black", pch = 21, asp=1, cex = .6)
plot(ann_new, add = TRUE, bg = "white", pch = 21, asp=1, cex = .6)
plot(G2, border = "darkred",add = TRUE, asp=1)
axis(1, at=seq(n3$lon-4, n3$lon+4, 2), labels=seq(0, 8, 2))
axis(1, at=seq(n3$lon-4, n3$lon+4, 0.5), labels=NA)
axis(2, at=seq(n3$lat-4, n3$lat+4, 2), labels=seq(0, 8, 2))
axis(2, at=seq(n3$lat-4, n3$lat+4, 0.5), labels=NA)
text(labels="a", cb$lon-7.5,cb$lat+7.5, adj=0, font = 2, cex = 1.5)
dev.off()
##### vioplots for news and olds ####
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

## 50x50 grid
r50 <- ggplot(data = tab05, aes(x = rec, y = sqrt(New), fill = rec)) +
  geom_violin( alpha = .3, scale =  "width", width = .6) +
  geom_boxplot(aes(y =sqrt(New), fill = rec), color = "black",
               width = 0.1, alpha = 0.4, guide = FALSE, outlier.shape = NA) +
  #expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis_d("site",labels = unique(tab2$rec)) +
  scale_fill_viridis_d("site",labels = unique(tab2$rec)) +
  theme_bw() +
  labs(x ="reef site", y = "Recruits per grid cell (root squared)")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text.x =  element_blank(),
        axis.ticks.x=  element_blank(),
        text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
r50
ggsave("figs/RecCount50_sqrt.png",width = 13, height = 9, units = 'cm')
a50 <- ggplot(data = tab05, aes(x = rec, y = sqrt(Old), fill = rec)) +
  geom_violin( alpha = .3, scale =  "width", width = .6) +
  geom_boxplot(color = "black",
               width = 0.1, alpha = 0.4, guide = FALSE, outlier.shape = NA) +
  #expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis_d("site",labels = unique(tab2$rec)) +
  scale_fill_viridis_d("site",labels = unique(tab2$rec)) +
  theme_bw() +
  labs(x ="reef site", y = "Adults per grid cell(square root)")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text.x =  element_blank(),
        axis.ticks.x=  element_blank(),
        text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
a50
ggsave("figs/AduCount50_sqrt.png",width = 13, height = 9, units = 'cm')

###### raw data pair plots ####

tab05$rec <- plyr::revalue(tab05$rec, c(CB19 = "CB", GT19 = "GT", SE19 = "SE", TB19 = "TB", TM19 = "TM", RS19 = "RS", L119 = "L1", N319  = "N3", HS19 = "HS"))
tab1$rec <- plyr::revalue(tab1$rec, c(CB19 = "CB", GT19 = "GT", SE19 = "SE", TB19 = "TB", TM19 = "TM", RS19 = "RS", L119 = "L1", N319  = "N3", HS19 = "HS"))
tab2$rec <- plyr::revalue(tab2$rec, c(CB19 = "CB", GT19 = "GT", SE19 = "SE", TB19 = "TB", TM19 = "TM", RS19 = "RS", L119 = "L1", N319  = "N3", HS19 = "HS"))

tab05pp<- tab05[,c("R","D","H","New", "Old", "mean_spats", "rec")]
tab1pp<- tab1[,c("R","D","H","New", "Old", "mean_spats", "rec")]
tab2pp<- tab2[,c("R","D","H","New", "Old", "mean_spats", "rec")]

pairplot05 <- ggpairs(data = tab05pp, 
        title = "Pair plot - 50cm grid",
        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=rec)), 
        diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=rec),
                    continuous = wrap("densityDiag", alpha=0.5)), 
        upper = list(combo = wrap("box_no_facet", alpha=0.5),
                     continuous = wrap("cor", size=4, alignPercent=0.8)),
        labeller = "label_parsed") + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
ggsave(file="figs/PairplotRaw50.png", plot=pairplot05, width=3, height=2,scale = 4)

pairplot1 <- ggpairs(data = tab1pp, 
        title = "Pair plot - 1m grid",
        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=rec)), 
        diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=rec),
                    continuous = wrap("densityDiag", alpha=0.5)), 
        upper = list(combo = wrap("box_no_facet", alpha=0.5),
                     continuous = wrap("cor", size=4, alignPercent=0.8)),
        labeller = "label_parsed") + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
ggsave(file="figs/PairplotRaw1m.png", plot=pairplot1, width=3, height=2,scale = 4)

pairplot2 <- ggpairs(data = tab2pp, 
        title = "Pair plot - 2m grid",
        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=rec)), 
        diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=rec),
                    continuous = wrap("densityDiag", alpha=0.5)), 
        upper = list(combo = wrap("box_no_facet", alpha=0.5),
                     continuous = wrap("cor", size=4, alignPercent=0.8)),
        labeller = "label_parsed") + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
ggsave(file="figs/PairplotRaw2m.png", plot=pairplot2, width=3, height=2,scale = 4)


## transform relevant data
tab05_trans <- as.data.frame(tab05 %>%
                                mutate(Olds = scale(Old),
                                       Oldl = log10(Old),
                                       Old_sqrt = sqrt(Old),
                                       Adus = scale(Adu),
                                       R_log10 = log10(R),
                                       Rs = as.numeric(scale(R)),
                                       Ds = as.numeric(scale(D)),
                                       mean_spatss = as.numeric(scale(mean_spats)),
                                       xCs = as.numeric(scale(xC)),
                                       yCs = as.numeric(scale(yC))))


tab1_trans <- as.data.frame(tab1 %>%
                               mutate(Olds = scale(Old),
                                      Oldl = log10(Old),
                                      Old_sqrt = sqrt(Old),
                                      Adus = scale(Adu),
                                      R_log10 = log10(R),
                                      Rs = as.numeric(scale(R)),
                                      Ds = as.numeric(scale(D)),
                                      mean_spatss = as.numeric(scale(mean_spats)),
                                      xCs = as.numeric(scale(xC)),
                                      yCs = as.numeric(scale(yC))))

tab2_trans <- as.data.frame(tab2 %>%
                               mutate(Olds = scale(Old),
                                      Oldl = log10(Old),
                                      Old_sqrt = sqrt(Old),
                                      Adus = scale(Adu),
                                      R_log10 = log10(R),
                                      Rs = as.numeric(scale(R)),
                                      Ds = as.numeric(scale(D)),
                                      mean_spatss = as.numeric(scale(mean_spats)),
                                      xCs = as.numeric(scale(xC)),
                                      yCs = as.numeric(scale(yC))))

### fit gam 

gam05 <- gam(data = tab05_trans, New ~ Old_sqrt + s(R_log10,D) + mean_spatss + s(rec, bs="re"), method="REML", family = quasipoisson())
summary(gam05)
tab_model(gam05)
gam.check(gam05, old.style = TRUE)
plot(residuals(gam05))
gam1 <- gam(data = tab1_trans, New ~ Old_sqrt + s(R_log10,D) + mean_spatss + s(rec, bs="re"), method="REML", family = poisson())
summary(gam1)
tab_model(gam1)
gam2 <- gam(data = tab2_trans, New ~ Old_sqrt + s(R_log10,D) + mean_spatss + s(rec, bs="re"), method="REML", family = poisson())
summary(gam2)
tab_model(gam2)

vis.gam(gam05, view=c("R_log10", "D"), type="response", plot.type="contour")
vis.gam(gam1, view=c("R_log10", "D"), type="response", plot.type="contour")
vis.gam(gam2, view=c("R_log10", "D"), type="response", plot.type="contour")


#### predictions ####
tab_model(gam05,gam1,gam2, show.df = TRUE,
          dv.labels = c("50x50cm grid", "1x1m grid","2x2m grid"),
          digits = 3)

grid.lines <- 100

D.pred.2 <- seq(min(tab2_trans$D), max(tab2_trans$D), length=grid.lines)
R.pred.2 <- seq(min(tab2_trans$R_log10), max(tab2_trans$R_log10), length=grid.lines)

mat <- matrix(1, length(R.pred.2), length(D.pred.2))
maty <- matrix(D.pred.2, length(R.pred.2), length(D.pred.2), byrow = TRUE)
matx <- matrix(R.pred.2, length(R.pred.2), length(D.pred.2))

mat.s <- expand.grid(R=R.pred.2, D=D.pred.2)
pred.df.2 <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Old_sqrt = rep(mean(tab2_trans$Old_sqrt), nrow(mat.s)),
                         mean_spatss = rep(mean(tab2_trans$mean_spatss), nrow(mat.s)),
                         xCs = rep(mean(tab2_trans$xCs), nrow(mat.s)),
                         yCs = rep(mean(tab2_trans$yCs), nrow(mat.s)),
                         rec = rep("NEW", nrow(mat.s)))

pred.df.2b <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Adu = rep(0, nrow(mat.s)),
                          mean_spats = rep(0, nrow(mat.s)),
                          rec = rep("NEW", nrow(mat.s)))

pred.s <- predict(gam2,pred.df.2)
pred.s <- predict(brm2_logsite,pred.df.2b, allow_new_levels=TRUE)
pred.df.2$pred <- pred.s
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.2 = xtabs(pred ~ R_log10 + D, data = pred.df.2)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))

pred.mat.store <- pred.mat.2  
pred.mat.2 <- pred.mat.store

pred.mat.2[!in.chull(matx, maty, 
                      tab2_trans$R_log10[chull(tab2_trans$R_log10, tab2_trans$D)], 
                      tab2_trans$D[chull(tab2_trans$R_log10, tab2_trans$D)])] <- NA


str(tab2_trans)

png("figs/pred2.png", width = 400, height = 400)
image(R.pred.2,D.pred.2, z = pred.mat.2, col=rev(hcl.colors(100, "reds",  alpha=0.7)), # axes=FALSE, 
      xlim=c(min(tab2_trans$R_log10), max(tab2_trans$R_log10)), 
      ylim=c(min(tab2_trans$D), max(tab2_trans$D)), xlab = expression("log"[10]*"(R)"), ylab = "D")

#axis(2, at=seq(1, 2.6, 0.2))
#axis(1, at=seq(0.75, 1.87, 0.2))
points(tab2_trans$R_log10, tab2_trans$D, col=rgb(0.5,0,0.3,0.1), pch=20)
contour(R.pred.2,D.pred.2, pred.mat.2, method = "edge", # levels=c(0, 0.1, 0.15, 0.2,0.25, 0.3, 0.4, 0.5, 0.7, 1, 1.5),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold"))
title("50x50cm grid")
dev.off()


D.pred.1 <- seq(min(tab1_trans$D), max(tab1_trans$D), length=grid.lines)
R.pred.1 <- seq(min(tab1_trans$R_log10), max(tab1_trans$R_log10), length=grid.lines)

mat <- matrix(1, length(R.pred.1), length(D.pred.1))
maty <- matrix(D.pred.1, length(R.pred.1), length(D.pred.1), byrow = TRUE)
matx <- matrix(R.pred.1, length(R.pred.1), length(D.pred.1))

mat.s <- expand.grid(R=R.pred.1, D=D.pred.1)
pred.df.1 <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Old_sqrt = rep(mean(tab1_trans$Old_sqrt), nrow(mat.s)),
                         mean_spatss = rep(mean(tab1_trans$mean_spatss), nrow(mat.s)),
                         xCs = rep(mean(tab1_trans$xCs), nrow(mat.s)),
                         yCs = rep(mean(tab1_trans$yCs), nrow(mat.s)),
                         rec = rep("NEW", nrow(mat.s)))

pred.df.1b <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Adu = rep(0, nrow(mat.s)),
                          mean_spats = rep(0, nrow(mat.s)),
                          rec = rep("NEW", nrow(mat.s)))

pred.s <- predict(gam1,pred.df.1)
pred.s <- predict(brm1_logsite,pred.df.1b, allow_new_levels=TRUE)
pred.df.1$pred <- pred.s
#write.csv(pred.df,"data/output/pred1.df.csv")

pred.mat.1 = xtabs(pred ~ R_log10 + D, data = pred.df.1)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))

pred.mat.store <- pred.mat.1  
pred.mat.1 <- pred.mat.store

pred.mat.1[!in.chull(matx, maty, 
                      tab1_trans$R_log10[chull(tab1_trans$R_log10, tab1_trans$D)], 
                      tab1_trans$D[chull(tab1_trans$R_log10, tab1_trans$D)])] <- NA


str(tab1_trans)

png("figs/pred1.png", width = 400, height = 400)
image(R.pred.1,D.pred.1, z = pred.mat.1, col=rev(hcl.colors(100, "reds",  alpha=0.7)), # axes=FALSE, 
      xlim=c(min(tab1_trans$R_log10), max(tab1_trans$R_log10)), 
      ylim=c(min(tab1_trans$D), max(tab1_trans$D)), xlab = expression("log"[10]*"(R)"), ylab = "D")

#axis(2, at=seq(1, 2.6, 0.2))
#axis(1, at=seq(0.75, 1.87, 0.2))
points(tab1_trans$R_log10, tab1_trans$D, col=rgb(0.5,0,0.3,0.1), pch=20)
contour(R.pred.1,D.pred.1, pred.mat.1, method = "edge", # levels=c(0, 0.1, 0.15, 0.2,0.25, 0.3, 0.4, 0.5, 0.7, 1, 1.5),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold"))
title("1x1cm grid")
dev.off()





D.pred.2 <- seq(min(tab2_trans$D), max(tab2_trans$D), length=grid.lines)
R.pred.2 <- seq(min(tab2_trans$R_log10), max(tab2_trans$R_log10), length=grid.lines)

mat <- matrix(1, length(R.pred.2), length(D.pred.2))
maty <- matrix(D.pred.2, length(R.pred.2), length(D.pred.2), byrow = TRUE)
matx <- matrix(R.pred.2, length(R.pred.2), length(D.pred.2))

mat.s <- expand.grid(R=R.pred.2, D=D.pred.2)
pred.df.2 <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Old_sqrt = rep(mean(tab2_trans$Old_sqrt), nrow(mat.s)),
                         mean_spatss = rep(mean(tab2_trans$mean_spatss), nrow(mat.s)),
                         xCs = rep(mean(tab2_trans$xCs), nrow(mat.s)),
                         yCs = rep(mean(tab2_trans$yCs), nrow(mat.s)),
                         rec = rep("NEW", nrow(mat.s)))

pred.df.2b <- data.frame(R_log10 = mat.s[,1], D = mat.s[,2], Adu = rep(0, nrow(mat.s)),
                          mean_spats = rep(0, nrow(mat.s)),
                          rec = rep("NEW", nrow(mat.s)))

pred.s <- predict(gam2,pred.df.2)
pred.s <- predict(brm2_logsite,pred.df.2b, allow_new_levels=TRUE)
pred.df.2$pred <- pred.s
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.2 = xtabs(pred ~ R_log10 + D, data = pred.df.2)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))

pred.mat.store <- pred.mat.2  
pred.mat.2 <- pred.mat.store

pred.mat.2[!in.chull(matx, maty, 
                      tab2_trans$R_log10[chull(tab2_trans$R_log10, tab2_trans$D)], 
                      tab2_trans$D[chull(tab2_trans$R_log10, tab2_trans$D)])] <- NA


str(tab2_trans)

png("figs/pred2.png", width = 400, height = 400)
image(R.pred.2,D.pred.2, z = pred.mat.2, col=rev(hcl.colors(100, "reds",  alpha=0.7)), # axes=FALSE, 
      xlim=c(min(tab2_trans$R_log10), max(tab2_trans$R_log10)), 
      ylim=c(min(tab2_trans$D), max(tab2_trans$D)), xlab = expression("log"[10]*"(R)"), ylab = "D")

#axis(2, at=seq(1, 2.6, 0.2))
#axis(1, at=seq(0.75, 1.87, 0.2))
points(tab2_trans$R_log10, tab2_trans$D, col=rgb(0.5,0,0.3,0.1), pch=20)
contour(R.pred.2,D.pred.2, pred.mat.2, method = "edge", # levels=c(0, 0.1, 0.15, 0.2,0.25, 0.3, 0.4, 0.5, 0.7, 1, 1.5),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold"))
title("50x50cm grid")
dev.off()D.pred.05 <- seq(min(tab05_trans$D), max(tab05_trans$D), length=grid.lines)
D.pred.05 <- seq(min(tab05_trans$D), max(tab05_trans$D), length=grid.lines)
D.pred.05 <- seq(min(tab05_trans$D), max(tab05_trans$D), length=grid.lines)
D.pred.05 <- seq(min(tab05_trans$D), max(tab05_trans$D), length=grid.lines)
D.pred.05 <- seq(min(tab05_trans$D), max(tab05_trans$D), length=grid.lines)
R.pred.05 <- seq(min(tab05_trans$R_log10), max(tab05_trans$R_log10), length=grid.lines)