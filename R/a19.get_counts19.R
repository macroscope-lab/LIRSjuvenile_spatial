### count annotations within gridded areas  

# upload library
library(terra)
library(tidyverse)
library(sf)
library(GISTools) #count.poly()
library(RColorBrewer)

## read annotations and central coordinates ####
ann<- st_read("data_raw/rr2019/annotations/Export_Output_all.shp")
ann[is.na(ann$Species),]$Species <-"Spp"
crs(ann)
coords <- read.csv("data_raw/rr2019/RR_gps_coords.csv")
c <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
c <- st_transform(c, crs(ann))
c_df <- cbind(as.data.frame(st_coordinates(c)),c) [,1:3]
names(c_df) <- c("x_rr","y_rr","site")
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

# prepare grid coordinates dataframes #####
grid1 <- data.frame() #will store bottom left and central coords of each 1mx1m square replicate within each RR
grid2 <- data.frame() #will store bottom left and central coords of each 2mx2m square replicate within each RR
grid4 <- data.frame() #will store bottom left and central coords of each 4mx4m square replicate within each RR

# populate grid1
L <- 1
gridding <- seq(0, 7, by = L)

for (site in c_df$site_full) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site_full == site]-4
  yb <- c_df$y_rr[c_df$site_full == site]-4
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
      write.csv(grid1, "data/grid1_coords19.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}

# populate grid2
L <- 2
gridding <- seq(0, 6, by = L)

for (site in c_df$site_full) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site_full == site]-4
  yb <- c_df$y_rr[c_df$site_full == site]-4
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
      write.csv(grid2, "data/grid2_coords19.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}

# populate grid4
L <- 4
gridding <- seq(0, 4, by = L)
for (site in c_df$site_full) {
  rep <- 1
  xb <- c_df$x_rr[c_df$site_full == site]-4
  yb <- c_df$y_rr[c_df$site_full == site]-4
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
      write.csv(grid4, "data/grid4_coords19.csv", row.names = FALSE)
      rep <- rep + 1
    }
  }
}


## prepare files for poly.count and count dots in each square of the grids ####

# split juveniles and adults
ann <- ann %>% filter(Status != "gone",
                      Species != "soft")

summary(as.factor(ann$Species))
sum(is.na(ann$Species))

ann_juv <- ann %>% filter(Status == "new" | Species == "Spp")
ann_adu <- ann %>% filter(Status != "new")


#grid1
c1 <- st_as_sf(grid1[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c1) <- st_crs(ann)

squares1 <- st_buffer(c1, dist = 0.5, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid1$n_adu <- poly.counts(ann_adu,squares1)
grid1$n_juv <- poly.counts(ann_juv,squares1)
write.csv(grid1, "data/grid1_counts19.csv", row.names = FALSE)


#grid2
c2 <- st_as_sf(grid2[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c2) <- st_crs(ann)

squares2 <- st_buffer(c2, dist = 1, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid2$n_adu <- poly.counts(ann_adu,squares2)
grid2$n_juv <- poly.counts(ann_juv,squares2)
write.csv(grid2, "data/grid2_counts19.csv", row.names = FALSE)


# grid4
c4 <- st_as_sf(grid4[,c("site","rep", "xc", "yc")],
               coords = c("xc", "yc"))
st_crs(c4) <- st_crs(ann)

squares4 <- st_buffer(c4, dist = 2, endCapStyle = "SQUARE") %>%
  st_cast("POLYGON")

grid4$n_adu <- poly.counts(ann_adu,squares4)
grid4$n_juv <- poly.counts(ann_juv,squares4)
write.csv(grid4, "data/grid4_counts19.csv", row.names = FALSE)


## plots ####

subset_sites <- c("lagoon_1",
                  "north_reef_3","gnarly_tree","corner_beach",
                  "horseshoe", "resort", "southeast",
                  "trimodal","turtle_beach")
for (rr in subset_sites) {
  filtered_squares1 <- squares1 %>% filter(site == rr)
  filtered_squares2 <- squares2 %>% filter(site == rr)
  filtered_squares4 <- squares4 %>% filter(site == rr)
  
  p <- ggplot() +
    geom_sf(data = ann_adu, color = "black") +
    geom_sf(data = ann_juv, aes(color = Observer)) +
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
