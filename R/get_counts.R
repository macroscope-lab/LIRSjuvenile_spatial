
library(raster)
library(tidyverse)
library(sf)
#library(rgeos)  
#library(rgdal) 
library(GISTools) #count.poly()
library(RColorBrewer)

## read coords and make grid files

coords <- read.csv("data_raw/sites.csv")
c2 <- st_as_sf(coords,
               coords = c("x", "y"),
               crs = 4326)
c <- st_transform(c2, CRS("+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
coords <- st_coordinates(c)
coord_t <- as.data.frame(cbind(site =c$Site, x = as.numeric(coords[,1]), y = as.numeric (coords[,2])))
coord_t$x <- as.numeric(coord_t$x)
coord_t$y <- as.numeric(coord_t$y)
# get bottom left coords of main square
coord_t$xb <- coord_t$x-4
coord_t$yb <- coord_t$y-4

#initialise grid coordinates dataframes #####
grid1 <- data.frame() #will store bottom left and central coords of each 1mx1m square replicate within each RR
grid2 <- data.frame() #will store bottom left and central coords of each 2mx2m square replicate within each RR
grid4 <- data.frame() #will store bottom left and central coords of each 4mx4m square replicate within each RR

# populate grid1
L <- 1
gridding <- seq(0, 7, by = L)

for (site in coord_t$site) {
  rep <- 1
  xb <- coord_t$xb[coord_t$site == site]
  yb <- coord_t$yb[coord_t$site == site]
  for (i in gridding) {
    for (j in gridding) {
      x0 <- xb + i
      y0 <- yb + j
      xc <- x0 + L / 2
      yc <- y0 + L / 2
      grid1 <-
        rbind(grid1,
              data.frame(
                site = site,
                rep = rep,
                x0 = x0,
                y0 = y0,
                xc = xc,
                yc = yc
              ))
      write.csv(grid1, "data/grid1_coords.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}

# populate grid2
L <- 2
gridding <- seq(0, 6, by = L)

for (site in coord_t$site) {
  rep <- 1
  xb <- coord_t$xb[coord_t$site == site]
  yb <- coord_t$yb[coord_t$site == site]
  for (i in gridding) {
    for (j in gridding) {
      x0 <- xb + i
      y0 <- yb + j
      xc <- x0 + L / 2
      yc <- y0 + L / 2
      grid2 <-
        rbind(grid2,
              data.frame(
                site = site,
                rep = rep,
                x0 = x0,
                y0 = y0,
                xc = xc,
                yc = yc
              ))
      write.csv(grid2, "data/grid2_coords.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}

# populate grid4
L <- 4
gridding <- seq(0, 4, by = L)
for (site in coord_t$site) {
  rep <- 1
  xb <- coord_t$xb[coord_t$site == site]
  yb <- coord_t$yb[coord_t$site == site]
  for (i in gridding) {
    for (j in gridding) {
      x0 <- xb + i
      y0 <- yb + j
      xc <- x0 + L / 2
      yc <- y0 + L / 2
      grid4 <-
        rbind(grid4,
              data.frame(
                site = site,
                rep = rep,
                x0 = x0,
                y0 = y0,
                xc = xc,
                yc = yc
              ))
      write.csv(grid4, "data/grid4_coords.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}

## read annotations ####
ann <- st_read("data_raw/annotations/lizard_rr_annotations_2022_nov.shp")


## prepare files fo poly.count ####

c4 <- st_as_sf(grid4[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c4) <- st_crs(ann)

squares4 <- st_buffer(c4, dist = 2, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON") 


### append counts on grid ####
grid4$countAdults <- poly.counts(ann,squares4)
write.csv(grid4, "data/grid4_counts.csv", row.names = FALSE)


