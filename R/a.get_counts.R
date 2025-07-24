### count annotations within gridded areas  

# upload library
library(terra)
library(tidyverse)
library(sf)
library(GISTools) #count.poly()
library(RColorBrewer)

## read annotations and central coordinates ####
ann<- st_read("data_raw/ann_shifted_manual.shp")
coords <- st_read("data_raw/dem_coords_nov2022_repr.shp")
ann22 <- st_read("data_raw/annotations/2022_04/lizard_rr_annotations_2022_april.shp")
ann22 <- st_transform(ann22, st_crs(ann))
crs(coords)
c <- st_transform(coords, crs(ann))
c_df <- cbind(as.data.frame(st_coordinates(c)),c) [,1:3]
names(c_df) <- c("x_rr","y_rr","site")
c_df <- c_df %>%
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

# prepare grid coordinates dataframes #####
grid1 <- data.frame() #will store bottom left and central coords of each 1mx1m square replicate within each RR
grid2 <- data.frame() #will store bottom left and central coords of each 2mx2m square replicate within each RR
grid4 <- data.frame() #will store bottom left and central coords of each 4mx4m square replicate within each RR

# populate grid1
L <- 1
gridding <- seq(0, 7, by = L)

for (site in c_df$site) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site == site]-4
  yb <- c_df$y_rr[c_df$site == site]-4
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

for (site in c_df$site) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site == site]-4
  yb <- c_df$y_rr[c_df$site == site]-4
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
for (site in c_df$site) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site == site]-4
  yb <- c_df$y_rr[c_df$site == site]-4
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


## prepare files for poly.count and count dots in each square of the grids ####

# split juveniles and adults
ann <- ann %>% filter(Status != "gone",
                      Species != "soft")
ann_juv <- ann %>% filter(ColonyID > max(ann22$ColonyID))
ann_adu <- ann %>% filter(ColonyID <= max(ann22$ColonyID))
#grid1
c1 <- st_as_sf(grid1[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c1) <- st_crs(ann)

squares1 <- st_buffer(c1, dist = 0.5, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid1$n_adu <- poly.counts(ann_adu,squares1)
grid1$n_juv <- poly.counts(ann_juv,squares1)
write.csv(grid1, "data/grid1_counts.csv", row.names = FALSE)


#grid2
c2 <- st_as_sf(grid2[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c2) <- st_crs(ann)

squares2 <- st_buffer(c2, dist = 1, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid2$n_adu <- poly.counts(ann_adu,squares2)
grid2$n_juv <- poly.counts(ann_juv,squares2)
write.csv(grid2, "data/grid2_counts.csv", row.names = FALSE)


# grid4
c4 <- st_as_sf(grid4[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c4) <- st_crs(ann)

squares4 <- st_buffer(c4, dist = 2, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid4$n_adu <- poly.counts(ann_adu,squares4)
grid4$n_juv <- poly.counts(ann_juv,squares4)
write.csv(grid4, "data/grid4_counts.csv", row.names = FALSE)


## plots ####

subset_sites <- c("cooks_path","easter_point","lagoon_2",
                  "north_reef_1","gnarly_tree","corner_beach",
                  "horseshoe", "north_reef_3", "resort", "southeast",
                  "trimodal","turtle_beach")
for (rr in subset_sites) {
  filtered_squares1 <- squares1 %>% filter(site == rr)
  filtered_squares2 <- squares2 %>% filter(site == rr)
  filtered_squares4 <- squares4 %>% filter(site == rr)
  
  p <- ggplot() +
    geom_sf(data = ann_adu, color = "black") +
    geom_sf(data = ann_juv, color = "blue") +
    geom_sf(data = filtered_squares1,fill = NA,color = "gold",linewidth = 1) +
    geom_sf(data = filtered_squares2,fill = NA,color = "orange",linewidth = 1) +
    geom_sf(data = filtered_squares4,fill = NA,color = "red",linewidth = 1) +
    theme_minimal() +
    ylim (ext(filtered_squares4)$ymin - 3, ext(filtered_squares4)$ymax + 3) +
    xlim(ext(filtered_squares4)$xmin - 3, ext(filtered_squares4)$xmax + 3) +
    ggtitle(rr)
  
  ggsave(paste0("figures/grids/",rr,"-ann_grid.png"), 
         plot = p, width = 6, height = 7, dpi = 300)
}


##### OLD #####
## read annotations and get site central coordinates ####
ann<- st_read("data_raw/annotations/lizard_rr_annotations_2022_nov.shp")
ann <- ann %>% filter(ann$Species != "soft")
max(ann$ColonyID) #42991

# some points in Vickis were labelled lagoon_1, which offset finding the centers.
# so taking L1 from previous annotation file

ann22 <- st_read("data_raw/annotations/2022_04/lizard_rr_annotations_2022_april.shp")
ann22 <- st_transform(ann22, st_crs(ann))
max(ann22$ColonyID) # 35179
ann22 <- ann22 %>% filter(Site != "lagoon_1")
ann_l1 <- ann %>% filter(Site == "lagoon_1")
ann2 <-rbind(ann22,ann_l1) 

# extract coordinates with features
coords <- st_coordinates(ann2)
points_df <- st_drop_geometry(ann2)
coords_df <- cbind(points_df, coords)

# Group by site and compute central x and y
centroids_df <- coords_df %>%
  group_by(Site) %>%
  summarize(
    xc = min(X) + (max(X)-min(X))/2,
    yc = min(Y) + (max(Y)-min(Y))/2,
    .groups = "drop"
  )

centroids_sf <- st_as_sf(centroids_df, coords = c("xc", "yc"), crs = st_crs(ann))


# prep for grid files

coords <- st_coordinates(centroids_sf)
coord_t <- as.data.frame(cbind(site =centroids_sf$Site, x = as.numeric(coords[,1]), y = as.numeric (coords[,2])))
coord_t$x <- as.numeric(coord_t$x)
coord_t$y <- as.numeric(coord_t$y)
# get bottom left coords of main square
coord_t$xb <- coord_t$x-4
coord_t$yb <- coord_t$y-4

# prepare grid coordinates dataframes #####
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


## prepare files for poly.count and count dots in each square of the grids ####

# split juveniles and adults
ann_juv <- ann %>% filter(ColonyID > max(ann22$ColonyID))
ann_adu <- ann %>% filter(ColonyID <= max(ann22$ColonyID),
                          Status != "gone")
#grid1
c1 <- st_as_sf(grid1[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c1) <- st_crs(ann)

squares1 <- st_buffer(c1, dist = 0.5, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid1$n_adu <- poly.counts(ann_adu,squares1)
grid1$n_juv <- poly.counts(ann_juv,squares1)
write.csv(grid1, "data/grid1_counts.csv", row.names = FALSE)


#grid2
c2 <- st_as_sf(grid2[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c2) <- st_crs(ann)

squares2 <- st_buffer(c2, dist = 1, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid2$n_adu <- poly.counts(ann_adu,squares2)
grid2$n_juv <- poly.counts(ann_juv,squares2)
write.csv(grid2, "data/grid2_counts.csv", row.names = FALSE)


# grid4
c4 <- st_as_sf(grid4[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c4) <- st_crs(ann)

squares4 <- st_buffer(c4, dist = 2, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid4$n_adu <- poly.counts(ann_adu,squares4)
grid4$n_juv <- poly.counts(ann_juv,squares4)
write.csv(grid4, "data/grid4_counts.csv", row.names = FALSE)


## plots ####

for (rr in unique(ann$Site)) {
  filtered_squares1 <- squares1 %>% filter(site == rr)
  filtered_squares2 <- squares2 %>% filter(site == rr)
  filtered_squares4 <- squares4 %>% filter(site == rr)
  
  p <- ggplot() +
    geom_sf(data = ann, color = "black") +
    geom_sf(data = ann_juv, color = "blue") +
    geom_sf(data = filtered_squares1,fill = NA,color = "gold",linewidth = 1) +
    geom_sf(data = filtered_squares2,fill = NA,color = "orange",linewidth = 1) +
    geom_sf(data = filtered_squares4,fill = NA,color = "red",linewidth = 1) +
    theme_minimal() +
    ylim (ext(filtered_squares4)$ymin - 3, ext(filtered_squares4)$ymax + 3) +
    xlim(ext(filtered_squares4)$xmin - 3, ext(filtered_squares4)$xmax + 3) +
    ggtitle(rr)
  
  ggsave(paste0("figures/grids/",rr,"-ann_grid.png"), 
         plot = p, width = 6, height = 7, dpi = 300)
}

