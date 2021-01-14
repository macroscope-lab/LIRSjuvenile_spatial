
library(raster)
library(tidyverse)
library(sp)
library(rgeos)  
library(rgdal) 
library(GISTools) #count.poly()
library(RColorBrewer)

## read coords and make grids with buffer

tab <- read.csv("other_output/grid2_theory_extra.csv", h = T)

summary(tab)
tab$xC <- tab$x0 +1
tab$yC <- tab$y0 +1

Ccoord2 <- tab[,c("xC", "yC")]
c2 <- SpatialPoints(Ccoord2)
tab2<-tab

tab1 <- read.csv("other_output/grid1_theory_extra.csv", h = T)

summary(tab1)
tab1$xC <- tab1$x0 +0.5
tab1$yC <- tab1$y0 +0.5

Ccoord1 <- tab1[,c("xC", "yC")]
c1 <- SpatialPoints(Ccoord1)



tab05 <- read.csv("other_output/grid05_theory_extra.csv", h = T)

summary(tab05)
tab05$xC <- tab05$x0 +0.25
tab05$yC <- tab05$y0 +0.25

Ccoord05 <- tab05[,c("xC", "yC")]
c05 <- SpatialPoints(Ccoord05)



## read annotations
ann <- readOGR("data_raw/annotations/Export_Output_all.shp")
ann_rec <- subset(ann, (grepl("rec", ann$Species) | grepl("R", ann$Species)) & ann$Species!= "Lrec")
unique(ann_rec$Species)

ann_adults <- subset(ann, (!grepl("rec", ann$Species) & !grepl("R", ann$Species))|ann$Species == "Lrec")
unique(ann_adults$Species)

ann_new <- subset(ann, ann$Status == "new")
ann_old <- subset(ann, ann$Status == "alive")

# read site rasters 
TM <- raster("data_raw/records/TM19.tif")
#TM <- projectRaster(TM, crs = crs(ref))
TM <- projectRaster(TM, 
                    crs = "+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
TB <- raster("data_raw/records/TB19.tif")
#TB <- projectRaster(TB, crs = crs(ref))
TB <- projectRaster(TB, 
                    crs = "+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
N3 <- raster("data_raw/records/N319.tif")
#N3 <- projectRaster(N3, crs = crs(ref))
N3 <- projectRaster(N3, 
                    crs = "+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
GT <- raster("data_raw/records/GT19.tif")
#GT <- projectRaster(GT, crs = crs(ref))
GT <- projectRaster(GT, 
                    crs = "+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
brewer.pal(n = 8, name = "BrBG")
# palette
ncols <- 100
newcol <- colorRampPalette(c("#8C510A", "#BF812D",
                             "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"))
cols <- newcol(ncols)  #apply the function to get 100 colours

# plot dem and rectangles
plot(TM, asp=1, col=cols, main = "Trimodal 2019")

# plot annotation
plot(ann_old, add = TRUE, col = "white")
plot(ann_new, add = TRUE, col = "black")

# plot grids
g2 <- list()
for(i in seq_along(c2)) {
  g2[i] <- gBuffer(c2[i], width=1, quadsegs=1, capStyle="SQUARE")
}
G2 <- do.call(bind, g2)
crs(G2) <- crs(ann)

# plot grids with buffer
gb2 <- list()
for(i in seq_along(c2)) {
  gb2[i] <- gBuffer(c2[i], width=1.10, quadsegs=1, capStyle="SQUARE")
}
GB2 <- do.call(bind, gb2)
crs(GB2) <- crs(ann)
plot(GB2, border = "white",add = TRUE)
plot(G2, border = "black",add = TRUE)


plot(N3, asp=1, col=cols, main = "North 3 2019")
plot(ann_adults, add = TRUE, col = "white")
plot(ann_rec, add = TRUE, col = "black")
plot(GB2, border = "white",add = TRUE)
plot(G2, border = "black",add = TRUE)

plot(GT, asp=1, col=cols, main = "Gnarly tree 2019")
plot(ann_adults, add = TRUE, col = "white")
plot(ann_rec, add = TRUE, col = "black")
plot(GB2, border = "white",add = TRUE)
plot(G2, border = "black",add = TRUE)

plot(TM, asp=1, col=cols, main = "Trimodal 2019")
plot(ann_adults, add = TRUE, col = "white")
plot(ann_rec, add = TRUE, col = "black")
plot(GB2, border = "white",add = TRUE)
plot(G2, border = "black",add = TRUE)

# count others in grid with buffer
countAdults = poly.counts(ann_adults,GB2)
dfAdults = data.frame(countAdults)

# count recs in grids
countRecs = poly.counts(ann_rec,G2)
dfRecs = data.frame(countRecs)

tab2$Rec <- dfRecs$countRecs
tab2$Adu <- dfAdults$countAdults


countnew2 = poly.counts(ann_new,G2)
dfnew2 = data.frame(countnew2)

countold2 = poly.counts(ann_old,GB2)
dfold2 = data.frame(countold2)

tab2$New <- dfnew2$countnew2
tab2$Old <- dfold2$countold2

# plot grids
g1 <- list()
for(i in seq_along(c1)) {
  g1[i] <- gBuffer(c1[i], width=0.5, quadsegs=1, capStyle="SQUARE")
}
G1 <- do.call(bind, g1)
crs(G1) <- crs(ann)

# plot grids with buffer
gb1 <- list()
for(i in seq_along(c1)) {
  gb1[i] <- gBuffer(c1[i], width=0.6, quadsegs=1, capStyle="SQUARE")
}
GB1 <- do.call(bind, gb1)
crs(GB1) <- crs(ann)

countRecs1 = poly.counts(ann_rec,G1)
dfRecs1 = data.frame(countRecs1)

countAdults1 = poly.counts(ann_adults,GB1)
dfAdults1 = data.frame(countAdults1)

tab1$Rec <- dfRecs1$countRecs1
tab1$Adu <- dfAdults1$countAdults1

countnew1 = poly.counts(ann_new,G1)
dfnew1 = data.frame(countnew1)

countold1 = poly.counts(ann_old,GB1)
dfold1 = data.frame(countold1)

tab1$New <- dfnew1$countnew1
tab1$Old <- dfold1$countold1



# plot grids
g05 <- list()
for(i in seq_along(c05)) {
  g05[i] <- gBuffer(c05[i], width=0.25, quadsegs=1, capStyle="SQUARE")
}
G05 <- do.call(bind, g05)
crs(G05) <- crs(ann)

# plot grids with buffer
gb05 <- list()
for(i in seq_along(c05)) {
  gb05[i] <- gBuffer(c05[i], width=0.26, quadsegs=1, capStyle="SQUARE")
}
GB05 <- do.call(bind, gb05)
crs(GB05) <- crs(ann)


countRecs05 = poly.counts(ann_rec,G05)
dfRecs05 = data.frame(countRecs05)

countAdults05 = poly.counts(ann_adults,GB05)
dfAdults05 = data.frame(countAdults05)

tab05$Recs <- dfRecs05$countRecs05
tab05$Adu <- dfAdults05$countAdults05

countnew05 = poly.counts(ann_new,G05)
dfnew05 = data.frame(countnew05)

countold05 = poly.counts(ann_old,GB05)
dfold05 = data.frame(countold05)

tab05$New <- dfnew05$countnew05
tab05$Old <- dfold05$countold05

#### add spats #####

spats <- read.csv("data_raw/tiles/settlement_2019-20_LIRS_VB.csv", head = TRUE)
spats <- as.data.frame(spats %>%
  group_by(site) %>%
  mutate(spats = sum(total), Aspat = sum(acroporidae)))
spats <- as.data.frame(spats %>%
                         group_by(site,tile_number) %>%
                         mutate(tile_spats = sum(total), tile_Aspat = sum(acroporidae)))
spats <- as.data.frame(spats %>%
                         distinct(site_abbreviation, tile_number,spats, Aspat, tile_spats, tile_Aspat))
spats <- as.data.frame(spats %>%
                         group_by(site_abbreviation) %>%
                         mutate(mean_spats = mean(tile_spats), mean_Aspat = mean(tile_Aspat)))
spats <- as.data.frame(spats %>%
                         distinct(site_abbreviation,mean_spats, mean_Aspat,spats, Aspat))
names(spats) [1]<- "rec"

# merge spats in 
tab05 <- merge(tab05, spats)
tab1 <- merge(tab1, spats)
tab2 <- merge(tab2, spats)

write.csv(tab05,"other_output/tab05_extra.csv")
write.csv(tab1,"other_output/tab1_extra.csv")
write.csv(tab2,"other_output/tab2_extra.csv")

