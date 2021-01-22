##get complexity ###
#  *modified from Torres-Pulliza et al.2020*  # 

## upload libraries  #install them first in case you don't have them already!
library(raster)
library(parallel)
library(rgdal)
library(sp)

## upload all this functions below

# Functions and packages

library(plot3D)
library(sgeostat)

# Load reference tif
ref <-  raster("data_raw/TB17.tif")

# Reef records

files <- dir("data_raw/records")
files <- unique(gsub("\\.tfw|\\.tif", "", files))
#files <- files[c(1,3,5,7,9,12,14,16,18)]
#
## run this loop that for each record gets heigth variation at different scales (necessary for computing D)
store05 <- data.frame()
store25 <- data.frame()
store1 <- data.frame()
store2 <- data.frame()
files

# read coords of the centers and transpose to right projection

coords <- read.csv("data_raw/RR_gps_coords.csv")
Ccoord <- coords[,c("x", "y")]
c2 <- SpatialPoints(Ccoord)
N3 <- raster("data_raw/records/N319.tif")
c2@proj4string <- crs(N3)
c <- spTransform(c2, CRS("+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
coords_tr<- as.data.frame(c)
Ccoords_tr<- cbind(coords$rec, coords$Site, coords_tr)
names(Ccoords_tr) <- c("rec","site","x","y")

for (rec in files) {
  # Load geotif for reef record
  data <- raster(paste0("data_raw/records/", rec, ".tif"))
  data <- projectRaster(data, crs = crs(ref))
  # Get lower corner of 8x8m bounding box; alternatively, prescribe these values.
  xb <- Ccoords_tr$x[which(Ccoords_tr$rec == rec)] - 4
  yb <- Ccoords_tr$y[which(Ccoords_tr$rec == rec)] - 4
  
  # Iterate through  0.5x0.5m quadrats (reps = 256) in reef record
  # Scope (extent), scales of variation, and resolution (grain)
  L <- 0.5 # Scope, 0.5 by 0.5 m reef patches
  scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude

  L0 <- min(scl) # Grain, resolution of processing ~ 7mm
  
  rep <- 1
  grid05 <- seq(0,7.5, by = 0.5)
  for (i in grid05) {
    for (j in grid05) {
      x0 <- xb + i
      y0 <- yb + j
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
      temp2 <- as(temp, 'SpatialGridDataFrame')
      R <- surfaceArea(temp2) / L^2
      D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
      store05 <- rbind(store05, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
      write.csv(store05, "other_output/grid05_notheory_centered.csv", row.names=FALSE)
      rep <- rep + 1
    }
  }
  
  # Iterate through  1x1 m quadrats (reps = 64) in reef record
  # Scope (extent), scales of variation, and resolution (grain)
  L <- 1 # Scope, 1 by 1 m reef patches
  scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  L0 <- min(scl) # Grain, resolution of processing ~ 16 mm
  
  rep <- 1
  grid1 <- seq(0,7, by = 1)
  for (i in grid1) {
    for (j in grid1) {
      x0 <- xb + i
      y0 <- yb + j
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
      temp2 <- as(temp, 'SpatialGridDataFrame')
      R <- surfaceArea(temp2) / L^2
      D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
      store1 <- rbind(store1, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
      write.csv(store1, "other_output/grid1_notheory_centered.csv", row.names=FALSE)
      rep <- rep + 1
    }
  }
  
  # Iterate through  2x2 m quadrats (reps = 16) in reef record
  # Scope (extent), scales of variation, and resolution (grain)
  L <- 2 # Scope, 2 by 2 m reef patches
  scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  L0 <- min(scl) # Grain, resolution of processing ~ 31 mm
  
  rep <- 1
  grid2 <- seq(0,6, by = 2)
  for (i in grid2) {
    for (j in grid2) {
      x0 <- xb + i
      y0 <- yb + j
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
      temp2 <- as(temp, 'SpatialGridDataFrame')
      R <- surfaceArea(temp2) / L^2
      D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
      store2 <- rbind(store2, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
      write.csv(store2, "other_output/grid2_notheory_centered.csv", row.names=FALSE)
      rep <- rep + 1
    }
  }

  L <- 0.25 # Scope, 0.5 by 0.5 m reef patches
  scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  L0 <- min(scl) # Grain, resolution of processing ~ 7mm
  
  rep <- 1
  grid25 <- seq(0,7.75, by = 0.25)
  for (i in grid25) {
    for (j in grid25) {
      x0 <- xb + i
      y0 <- yb + j
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
      temp2 <- as(temp, 'SpatialGridDataFrame')
      R <- surfaceArea(temp2) / L^2
      D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
      store25 <- rbind(store25, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
      write.csv(store25, "other_output/grid25_notheory_centered.csv", row.names=FALSE)
      rep <- rep + 1
    }
  }
}

files

#### extra
store05 <- data.frame()
store1 <- data.frame()
store2 <- data.frame()
files

for (rec in files) {
  # Load geotif for reef record
  data <- raster(paste0("data_raw/records/", rec, ".tif"))
  data <- projectRaster(data, crs = crs(ref))
  # Get lower corner of 8x8m bounding box; alternatively, prescribe these values.
  xb <- Ccoords_tr$x[which(Ccoords_tr$rec == rec)] - 4
  yb <- Ccoords_tr$y[which(Ccoords_tr$rec == rec)] - 4
  
  x1 <- xb - 2
  y1 <- yb + 2
  x2 <- xb - 2
  y2 <- yb + 4
  x3 <- xb + 2
  y3 <- yb + 8
  x4 <- xb + 4
  y4 <- yb + 8
  x5 <- xb + 8
  y5 <- yb + 4
  x6 <- xb + 8
  y6 <- yb + 2
  x7 <- xb + 4
  y7 <- yb - 2
  x8 <- xb + 2
  y8 <- yb - 2
  
  x <- c(x1, x2, x3, x4, x5, x6, x7, x8)
  y <- c(y1, y2, y3, y4, y5, y6, y7, y8)
  
  for (n in 1:8) {
    xb <- x[n]
    yb <- y[n]
    
    # Iterate through  0.5x0.5m quadrats (reps = 256) in reef record
    # Scope (extent), scales of variation, and resolution (grain)
    L <- 0.5 # Scope, 0.5 by 0.5 m reef patches
    scl <-
      L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
    L0 <- min(scl) # Grain, resolution of processing ~ 7mm
    
    rep <- 257
    grid05 <- seq(0, 1.5, by = 0.5)
    for (i in grid05) {
      for (j in grid05) {
        x0 <- xb + i
        y0 <- yb + j
        temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
        H <-
          max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
        temp2 <- as(temp, 'SpatialGridDataFrame')
        R <- surfaceArea(temp2) / L ^ 2
        D <-
          3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0)
        store05 <-
          rbind(store05,
                data.frame(
                  rec = rec,
                  rep = rep,
                  x0 = x0,
                  y0 = y0,
                  H = H,
                  R = R,
                  D = D
                ))
        write.csv(store05,
                  "other_output/grid05_theory_extra_centered.csv",
                  row.names = FALSE)
        rep <- rep + 1
      }
    }
    
    # Iterate through  1x1 m quadrats (reps = 64) in reef record
    # Scope (extent), scales of variation, and resolution (grain)
    L <- 1 # Scope, 1 by 1 m reef patches
    scl <-
      L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
    L0 <- min(scl) # Grain, resolution of processing ~ 16 mm
    
    rep <- 577
    grid1 <- seq(0, 1, by = 1)
    for (i in grid1) {
      for (j in grid1) {
        x0 <- xb + i
        y0 <- yb + j
        temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
        H <-
          max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
        temp2 <- as(temp, 'SpatialGridDataFrame')
        R <- surfaceArea(temp2) / L ^ 2
        D <-
          3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0)
        store1 <-
          rbind(store1,
                data.frame(
                  rec = rec,
                  rep = rep,
                  x0 = x0,
                  y0 = y0,
                  H = H,
                  R = R,
                  D = D
                ))
        write.csv(store1,
                  "other_output/grid1_theory_extra_centered.csv",
                  row.names = FALSE)
        rep <- rep + 1
      }
    }
    
    # Iterate through  2x2 m quadrats (reps = 16) in reef record
    # Scope (extent), scales of variation, and resolution (grain)
    L <- 2 # Scope, 2 by 2 m reef patches
    scl <-
      L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
    L0 <- min(scl) # Grain, resolution of processing ~ 31 mm
    
    rep <- 145
    grid2 <- seq(0, 1, by = 2)
    x0 <- xb
    y0 <- yb
    temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
    H <-
      max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
    temp2 <- as(temp, 'SpatialGridDataFrame')
    R <- surfaceArea(temp2) / L ^ 2
    D <-
      3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0)
    store2 <-
      rbind(store2,
            data.frame(
              rec = rec,
              rep = rep,
              x0 = x0,
              y0 = y0,
              H = H,
              R = R,
              D = D
            ))
    write.csv(store2, "other_output/grid2_theory_extra_centered.csv", row.names = FALSE)
    rep <- rep + 1
  }
  
  L <- 0.25 # Scope, 2 by 2 m reef patches
  scl <-
    L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  L0 <- min(scl) # Grain, resolution of processing ~ 31 mm
  
  rep <- 9217
  grid25 <- seq(0, 1.75, by = 0.25)
  for (i in grid1) {
    for (j in grid1) {
      x0 <- xb + i
      y0 <- yb + j
  temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
  H <-
    max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
  temp2 <- as(temp, 'SpatialGridDataFrame')
  R <- surfaceArea(temp2) / L ^ 2
  D <-
    3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0)
  store2 <-
    rbind(store2,
          data.frame(
            rec = rec,
            rep = rep,
            x0 = x0,
            y0 = y0,
            H = H,
            R = R,
            D = D
          ))
  write.csv(store2, "other_output/grid2_theory_extra_centered.csv", row.names = FALSE)
  rep <- rep + 1
  }

  }
}


## add extra quadrats

  
plot(store05$R~grid05$R_theory,xlab="R_theory", ylab="R",)
plot(store05$D~grid05$D_theory,xlab="D from R_theory", ylab="D from R",)
plot(store05$H~grid05$H)

plot(store1$R~tab1$R_theory,xlab="R_theory", ylab="R",)
plot(store1$D~tab1$D_theory,xlab="D from R_theory", ylab="D from R",)
plot(store1$H~tab1$H)

plot(store2$R~tab2$R_theory,xlab="R_theory", ylab="R",)
plot(store2$D~tab2$D_theory,xlab="D from R_theory", ylab="D from R",)
plot(store2$H~tab2$H)

## Fit more quadrats
