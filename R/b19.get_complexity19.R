##get complexity ###
# get complexity values for grids of different sizes

## upload libraries
library(terra)
library(parallel)
library(sf)
library(habtools)
library(raster)
library(dplyr)

# Load reference tif and centers
ref <- rast("data/dem_reprojected_resampled/TM19_repr.tif")
coords <- read.csv("data_raw/rr2019/RR_gps_coords.csv")
c <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
c <- st_transform(c, crs(ref))
c_df <- cbind(as.data.frame(st_coordinates(c)),c) [,1:3]
names(c_df) <- c("x_rr","y_rr","site")
subset_sites <- c(
  "lagoon_1",
  "north_reef_3",
  "gnarly_tree",
  "corner_beach",
  "horseshoe",
  "resort",
  "southeast",
  "trimodal",
  "turtle_beach"
)
c_df <- c_df %>%
  mutate(site_full = case_when(
    site == "CP" ~ "cooks_path",
    site == "CB19" ~ "corner_beach",
    site == "EP" ~ "easter_point",
    site == "GT19" ~ "gnarly_tree",
    site == "HS19" ~ "horseshoe",
    site == "L119" ~ "lagoon_1",
    site == "L219" ~ "lagoon_2",
    site == "LH" ~ "lizard_head",
    site == "NML" ~ "no_mans_land",
    site == "NP" ~ "north_of_paradise",
    site == "N1" ~ "north_reef_1",
    site == "N2" ~ "north_reef_2",
    site == "N319" ~ "north_reef_3",
    site == "OS" ~ "osprey",
    site == "RS19" ~ "resort",
    site == "SI" ~ "south_island",
    site == "SE19" ~ "southeast",
    site == "TM19" ~ "trimodal",
    site == "TB19" ~ "turtle_beach",
    site == "VK19" ~ "vickis",
    site == "WM" ~ "washing_machine")) %>%
  filter(site_full %in% subset_sites)

c <- st_as_sf(c_df, coords = c("x_rr", "y_rr"), crs = crs(ref))


## read file list
files <- dir("data_raw/rr2019/records")
files <- unique(sub("\\..*", "", files))
files

## resample and reproject ####
# resample and export dem at right lizard meter reprojection
res <- data.frame()
target_res <- 0.005

for (rr in files) {
  data <- rast(paste0("data_raw/rr2019/records/", rr, ".tif"))
  data <- terra::project(data, crs(ann))
  site <- rr
  resx <- res(data)[1]
  resy <- res(data) [2]
  ext_rr <- ext(data)
  target <- rast(
    extent = ext_rr,
    resolution = c(target_res, target_res),
    crs = crs(data)
  )
  rr_resampled <- resample(data, target, method = "bilinear")
  resx_re <- res(rr_resampled)[1]
  resy_re <- res(rr_resampled) [2]
  res <- rbind(res,cbind(site = site, resx = resx, resy = resy,
                         resx_re = resx_re, resy_re = resy_re))
  writeRaster(rr_resampled, paste0("data/dem_reprojected_resampled/",rr,"_repr.tif"))
}

write.csv(res,"outputs/raw_dem_reproj_resolutions2019.csv")

## get rdh ####
#seems fine - get tready to loop through
files <- dir("data/dem_reprojected_resampled")
files <- unique(sub("_.*", "", files))
files <- files[grepl("19", files)]

## run this loop for each scale
store50 <- data.frame()
store1 <- data.frame()
store2 <- data.frame()
store4 <- data.frame()

# loop through all rr
#rr <- "CP"
for (rr in files) {
  # Load geotif for reef record
  data <- rast(paste0("data/dem_reprojected_resampled/", rr, "_repr.tif"))
  center <- c[c$site==rr,]
  x0 <- as.numeric(st_coordinates(center)[,1])
  y0 <- as.numeric(st_coordinates(center)[,2])
  rrs <- dem_crop(raster(data),x0,y0,8)
  
  dem_list1 <- dem_split(rrs,habtools::extent(raster(rrs))/8)
  results1 <- lapply(dem_list1,rdh)
  result1_df <- do.call(rbind, results1)
  result1_df$site <- rr
  result1_df$rep <- c(1:64)
  store1 <- rbind(store1, result1_df)
  write.csv(store1, "data/grid1_rdh19.csv")
  
  dem_list2 <- dem_split(rrs,habtools::extent(raster(rrs))/4)
  results2 <- lapply(dem_list2,rdh)
  result2_df <- do.call(rbind, results2)
  result2_df$site <- rr
  result2_df$rep <- c(1:16)
  store2 <- rbind(store2, result2_df)
  write.csv(store2, "data/grid2_rdh19.csv")
  
  dem_list4 <- dem_split(rrs,habtools::extent(raster(rrs))/2)
  results4 <- lapply(dem_list4,rdh)
  result4_df <- do.call(rbind, results4)
  result4_df$site <- rr
  result4_df$rep <- c(1,3,2,4)
  store4 <- rbind(store4, result4_df)
  write.csv(store4, "data/grid4_rdh19.csv")
  
}

# fix rep ID to match rep ID in the count dataset

# S=1m
#store1 <- read.csv("data/grid1_rdh19.csv")
store1$rep_ht<- store1$rep

tr <- data.frame(old = seq(1:64))
tr$factor <- c(seq(0,7),seq(-1,6),seq(-2,5),seq(-3,4),
               seq(-4,3),seq(-5,2),seq(-6,1),seq(-7,0))
tr$new <- tr$old + (7*tr$factor)
store1$rep <- rep(tr$new,9)

write.csv(store1, "data/grid1_rdh19.csv")

# S=2m
#store2 <- read.csv("data/grid2_rdh19.csv")
store2$rep_ht<- store2$rep

tr <- data.frame(old = seq(1:16))
tr$factor <- c(seq(0,3),seq(-1,2),seq(-2,1),seq(-3,0))
tr$new <- tr$old + (3*tr$factor)
store2$rep <- rep(tr$new,9)

write.csv(store2, "data/grid2_rdh19.csv")

# S=4m
# coded different as reps were already fixed in the loop above
#store4 <- read.csv("data/grid4_rdh19.csv")

tr <- data.frame(old = c(1,3,2,4))
tr$new <- c(1:4)
store4$rep_ht<- rep(tr$new,9)

write.csv(store4, "data/grid4_rdh19.csv")
