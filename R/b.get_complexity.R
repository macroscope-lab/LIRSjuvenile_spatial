##get complexity ###
# get complexity values for grids of different sizes

## upload libraries
library(terra)
library(parallel)
library(sf)
library(habtools)
library(raster)

# Load reference tif and centers
rr <- rast("data_raw/records/OS.tif")
ann<- st_read("data_raw/annotations/lizard_rr_annotations_2022_nov.shp")
coords <- st_read("data_raw/dem_coords_nov2022.shp")
c_rep <- st_transform(coords, crs(ann))
# a plot
plot(rr)
points(c_rep)

## read file list
files <- dir("data_raw/records")
files <- unique(sub("\\..*", "", files))


## resample and reproject ####
# resample and export dem at right lizard meter reprojection
res <- data.frame()
target_res <- 0.005
for (rr in files) {
  data <- rast(paste0("data_raw/records/", rr, ".tif"))
  data <- project(data, crs(ann))
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
write.csv(res[res$site!="HS2",], "outputs/raw_dem_reproj_resolutions.csv")


## get rdh ####
#seems fine - get tready to loop through
files <- dir("data/dem_reprojected_resampled")
files <- unique(sub("_.*", "", files))

## run this loop for each scale
store1 <- data.frame()
store2 <- data.frame()
store4 <- data.frame()

## put centres in meter ref
c <- st_transform(coords, crs(ann))

# loop through all rr
#rr <- "CP"
for (rr in files) {
  # Load geotif for reef record
  data <- rast(paste0("data/dem_reprojected_resampled/", rr, "_repr.tif"))
  center <- c[c$rec==rr,]
  x0 <- as.numeric(st_coordinates(center)[,1])
  y0 <- as.numeric(st_coordinates(center)[,2])
  rrs <- dem_crop(raster(data),x0,y0,8)
  
  dem_list1 <- dem_split(rrs,habtools::extent(raster(rrs))/8)
  results1 <- lapply(dem_list1,rdh)
  result1_df <- do.call(rbind, results1)
  result1_df$site <- rr
  result1_df$rep <- c(1:64)
  store1 <- rbind(store1, result1_df)
  write.csv(store1, "data/grid1_rdh.csv")
  
  dem_list2 <- dem_split(rrs,habtools::extent(raster(rrs))/4)
  results2 <- lapply(dem_list2,rdh)
  result2_df <- do.call(rbind, results2)
  result2_df$site <- rr
  result2_df$rep <- c(1:16)
  store2 <- rbind(store2, result2_df)
  write.csv(store2, "data/grid2_rdh.csv")
  
  dem_list4 <- dem_split(rrs,habtools::extent(raster(rrs))/2)
  results4 <- lapply(dem_list4,rdh)
  result4_df <- do.call(rbind, results4)
  result4_df$site <- rr
  result4_df$rep <- c(1,3,2,4)
  store4 <- rbind(store4, result4_df)
  write.csv(store4, "data/grid4_rdh.csv")
  
}

# fix rep ID to match rep ID in the count dataset

store50$rep_ht<- store50$rep
tr <- data.frame(old = seq(1:256))
tr$factor <- c(seq(0,15),seq(-1,14),seq(-2,13),seq(-3,12),
seq(-4,11),seq(-5,10),seq(-6,9),seq(-7,8),
seq(-8,7),seq(-9,6),seq(-10,5),seq(-11,4),
seq(-12,3),seq(-13,2),seq(-14,1),seq(-15,0))
tr$new <- tr$old + (7*tr$factor)
tr$new <- tr$old + (15*tr$factor)
store50$rep <- rep(tr$new,9)
write.csv(store50, "data/grid50_rdh19.csv")


# S=1m
#store1 <- read.csv("data/grid1_rdh.csv")
store1$rep_ht<- store1$rep

tr <- data.frame(old = seq(1:64))
tr$factor <- c(seq(0,7),seq(-1,6),seq(-2,5),seq(-3,4),
               seq(-4,3),seq(-5,2),seq(-6,1),seq(-7,0))
tr$new <- tr$old + (7*tr$factor)
store1$rep <- rep(tr$new,18)

write.csv(store1, "data/grid1_rdh.csv")

# S=2m
#store2 <- read.csv("data/grid2_rdh.csv")
store2$rep_ht<- store2$rep

tr <- data.frame(old = seq(1:16))
tr$factor <- c(seq(0,3),seq(-1,2),seq(-2,1),seq(-3,0))
tr$new <- tr$old + (3*tr$factor)
store2$rep <- rep(tr$new,18)

write.csv(store2, "data/grid2_rdh.csv")

# S=4m
# coded different as reps were already fixed in the loop above
#store4 <- read.csv("data/grid4_rdh.csv")

tr <- data.frame(old = c(1,3,2,4))
tr$new <- c(1:4)
store4$rep_ht<- rep(tr$new,18)

write.csv(store4, "data/grid4_rdh.csv")

# compare site center coords with annotation file ####
c_rep <- st_transform(coords, crs(ann))
c_rr <- cbind(as.data.frame(st_coordinates(c_rep)),c_rep) [,1:3]
names(c_rr) <- c("x_rr","y_rr","site")
c_ann <- as.data.frame(unique(dt1[,c("site","xs","ys")]))
names(c_ann) <- c("site","x_ann","y_ann")
c_comp <- merge(c_ann,c_rr)
c_comp <- c_comp %>%
  mutate(vector_x = x_rr - x_ann,
         vector_y = y_rr - y_ann)
shift_df <- c_comp[,c("site", "vector_x", "vector_y")]

shift_df <- shift_df %>%
  mutate(site = case_when(
    site == "CP" ~ "cooks_path",
    site == "CB" ~ "corner_beach",
    site == "EP" ~ "easter_point",
    site == "GT" ~ "gnarly_tree",
    site == "HS" ~ "horseshoe",
    site == "L1" ~ "lagoon_1",
    site == "L2" ~ "lagoon_2",
    site == "LH" ~ "lizard_head",
    site == "NM" ~ "no_mans_land",
    site == "NP" ~ "north_of_paradise",
    site == "N1" ~ "north_reef_1",
    site == "N2" ~ "north_reef_2",
    site == "N3" ~ "north_reef_3",
    site == "OS" ~ "osprey",
    site == "RS" ~ "resort",
    site == "SI" ~ "south_island",
    site == "SE" ~ "southeast",
    site == "TM" ~ "trimodal",
    site == "TB" ~ "turtle_beach",
    site == "VK" ~ "vickis",
    site == "WM" ~ "washing_machine",
  ))

names(shift_df) <- c("Site", "vector_x", "vector_y")

ann_shifted <- ann[!ann$Site %in% c("north_reef_2","washing_machine","lizard_head", "no_mans_land"),] %>%
  left_join(shift_df, by = "Site")

coords <- st_coordinates(ann_shifted)

# Apply the shift
coords_shifted <- coords
coords_shifted[,1] <- coords[,1] + ann_shifted$vector_x
coords_shifted[,2] <- coords[,2] + ann_shifted$vector_y

ann_shifted$geometry <- st_sfc(lapply(1:nrow(coords_shifted), function(i) {
  st_point(coords_shifted[i, ])
}), crs = st_crs(ann_shifted))


## loop for plotting translated points
for (rr in files) {
  # Load geotif for reef record
  data <- rast(paste0("data_raw/records/", rr, ".tif"))
  data <- project(data, crs(ann))
  rr_df <- as.data.frame(data, xy = TRUE, na.rm = TRUE)
  raster_var <- names(data)[1]
  ext <- ext(data)
  # Step 2: Plot using ggplot
  p <- ggplot() +
    geom_raster(data = rr_df, aes(x = x, y = y, fill = !!sym(raster_var))) +
    geom_sf(data = ann_shifted, color = "red", size = 1) +
    geom_sf(data = ann, color = "blue", size = 1) +
    coord_sf(
      xlim = c(ext[1], ext[2]),  # xmin, xmax
      ylim = c(ext[3], ext[4]),  # ymin, ymax
      expand = FALSE
    ) +
    theme_minimal() +
    labs(title = paste0(rr, " raster with Shifted Points in red"), fill = raster_var)
  
  ggsave(paste0("figures/grids/",rr,"-ann_shifted.png"), 
         plot = p, width = 10, height = 10, dpi = 400)
}


st_write(ann_shifted, "data_raw/ann_shifted.shp", overwrite = TRUE)



test <- dir("data_raw/orthos")
test <- unique(sub("\\..*", "", test))
for (r in test) {
  # Load geotif for reef record
  data <- rast(paste0("data_raw/orthos/", r, ".tif"))
  data <- project(data, crs(ann))
  rr_df <- as.data.frame(data, xy = TRUE, na.rm = TRUE)
  raster_var <- names(data)[1]
  ext <- ext(data)
  # Step 2: Plot using ggplot
  # p <- ggplot() +
  #   geom_raster(data = rr_df, aes(x = x, y = y, fill = !!sym(raster_var))) +
  #   geom_sf(data = ann_shifted, color = "red", size = 1) +
  #   geom_sf(data = ann, color = "blue", size = 1) +
  #   coord_sf(
  #     xlim = c(ext[1], ext[2]),  # xmin, xmax
  #     ylim = c(ext[3], ext[4]),  # ymin, ymax
  #     expand = FALSE
  #   ) +
  #   theme_minimal() +
  #   labs(title = paste0(rr, " raster with Shifted Points in red"), fill = raster_var)
  # 
  # ggsave(paste0("figures/grids/",r,"-orthoproj_shifted.png"), 
  #        plot = p, width = 10, height = 10, dpi = 400)
  writeRaster(data, paste0("data_raw/orthos-repr/",r,"_repr.tif"), overwrite = TRUE)
}

### check ####
#rdh4 <- read.csv("data/grid4_rdh.csv")
rdh1 <- store1

names(rdh1)
summary(rdh1)
rdh1[c(tail(order(rdh1$R))),]
rdh1[c(head(order(rdh1$R))),]

rdh1[c(tail(order(rdh1$D))),]
rdh1[c(head(order(rdh1$D))),]

rdh1[c(tail(order(rdh1$H))),]
rdh1[c(head(order(rdh1$H))),]

rdh2[c(tail(order(rdh2$R))),]
rdh2[c(head(order(rdh2$R))),]

rdh2[c(tail(order(rdh2$D))),]
rdh2[c(head(order(rdh2$D))),]

rdh2[c(tail(order(rdh2$H))),]
rdh2[c(head(order(rdh2$H))),]

rdh4[c(tail(order(rdh4$R))),]
rdh4[c(head(order(rdh4$R))),]

rdh4[c(tail(order(rdh4$D))),]
rdh4[c(head(order(rdh4$D))),]

rdh4[c(tail(order(rdh4$H))),]
rdh4[c(head(order(rdh4$H))),]

 