##get complexity ###
# get complexity values for grids of different sizes

## upload libraries
library(terra)
library(parallel)
library(sf)
library(habtools)

# Load reference tif and centers
rr <- rast("data_raw/records/SI.tif")
coords <- st_read("data_raw/dem_coords_nov2022.shp")
c_rep <- st_transform(coords, crs(rr))
ann<- st_read("data_raw/annotations/lizard_rr_annotations_2022_nov.shp")
# a plot
plot(rr)
points(c_rep)

#seems fine - get tready to loop through
## read file list
files <- dir("data_raw/records")
files <- unique(sub("\\..*", "", files))

## run this loop for each scale
store1 <- data.frame()
store2 <- data.frame()
store4 <- data.frame()

## put centres in meter ref
c <- st_transform(coords, crs(ann))

## a plot to test projection change will work
rrp <- project(rr, crs(ann))
plot(rrp)
points(c)

# loop through all rr
#rr <- "CP"
for (rr in files) {
  # Load geotif for reef record
  data <- rast(paste0("data_raw/records/", rr, ".tif"))
  data <- project(data, crs(ann))
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

#### some checks
# resolutions
res <- data.frame()
for (rr in files) {
  data <- rast(paste0("data_raw/records/", rr, ".tif"))
  site <- rr
  resx <- res(data)[1]
  resy <- res(data) [2]
  res <- rbind(res,cbind(site = site, resx = resx, resy = resy))
}
write.csv(res[res$site!="HS2",], "outputs/raw_dem_resolutions.csv")

# compare site center coords ####
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

writeOGR
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

st_write(ann_shifted, "data_raw/ann_shifted.shp")
#it looks OK
