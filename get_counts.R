### get counts ####
# upload libraries
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

### functions ####
create_grid <- function(c_df, L, output_file = "outputs/gridX_coords.csv") {
  gridding <- seq(0, 7.5, by = L)
  grid <- c_df %>%
    select(site, x_rr, y_rr) %>%
    distinct() %>%
    # create all combinations of i and j for each site
    crossing(i = gridding, j = gridding) %>%
    mutate(
      xb = x_rr - 4,
      yb = y_rr - 4,
      x0 = xb + i,
      y0 = yb + j,
      xc = x0 + L / 2,
      yc = y0 + L / 2
    ) %>%
    group_by(site) %>%
    mutate(rep = row_number()) %>%
    ungroup() %>%
    select(site, rep, x0, y0, xc, yc)
  write.csv(grid, output_file, row.names = FALSE)
  return(grid)
}
grid_counts <- function(grid_list,
                        ann_adu,
                        ann_juv,
                        square_sizes,
                        output_file = "outputs/grid_countsXX.csv") {
  # function to process each grid
  process_one <- function(grid, square_size) {
    pts <- grid %>%
      select(site, rep, xc, yc) %>%
      st_as_sf(coords = c("xc", "yc"), remove = FALSE)
    st_crs(pts) <- st_crs(ann)
    half_size <- square_size / 2
    squares <- pts %>%
      st_buffer(dist = half_size, endCapStyle = "SQUARE") %>%
      st_cast("POLYGON")
    # counts
    n_adu <- lengths(st_intersects(squares, ann_adu))
    n_juv <- lengths(st_intersects(squares, ann_juv))
    squares <- squares %>%
      mutate(n_adu = n_adu,
             n_juv = n_juv,
             scale = square_size)
    return(squares)
  }
  # apply across all grids
  all_squares <- map2_dfr(grid_list, square_sizes, process_one)
  #save non-geometry version if needed
  out_df <- all_squares %>%
    st_drop_geometry()
  write.csv(out_df, output_file, row.names = FALSE)
  return(all_squares)
}
plot_grids_by_site <- function(all_squares,
                               ann_adu,
                               ann_juv,
                               buffer = 3 ,
                               path_year = 0) {
  for (rr in unique(all_squares$site)) {
    sq <- all_squares %>% filter(site == rr)
    bb <- st_bbox(sq %>% filter(scale == max(scale)))
    p <- ggplot() +
      geom_sf(data = ann_adu, color = "black") +
      geom_sf(data = ann_juv, color = "blue") +
      geom_sf(
        data = sq %>% filter(scale == 0.5),
        fill = NA,
        color = "green",
        linewidth = 0.8
      ) +
      geom_sf(
        data = sq %>% filter(scale == 1),
        fill = NA,
        color = "gold",
        linewidth = 1
      ) +
      geom_sf(
        data = sq %>% filter(scale == 2),
        fill = NA,
        color = "orange",
        linewidth = 1
      ) +
      geom_sf(
        data = sq %>% filter(scale == 4),
        fill = NA,
        color = "red",
        linewidth = 1
      ) +
      coord_sf(
        xlim = c(bb["xmin"] - buffer, bb["xmax"] + buffer),
        ylim = c(bb["ymin"] - buffer, bb["ymax"] + buffer)
      ) +
      ggtitle(rr) +
      theme_minimal()
    ggsave(
      paste0("figures/grids/", rr, "_grid_", path_year, ".png"),
      plot = p,
      width = 6,
      height = 7,
      dpi = 300
    )
  }
}
### get counts per each year #####

# everything needs to be transformed in meters projections, same as done for maps
ref <- rast("data/dem_reprojected_resampled/TM19_repr.tif")
ann <- st_read("data/rr_annotations_lirs_2019.shp")
ann[is.na(ann$Species), ]$Species <- "Spp"
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
crs(ann)
ann <- st_transform(ann, crs(ref))
coords <- read.csv("data/rr19_gps_coords.csv")
c <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
c <- st_transform(c, crs(ref))
c_df <- cbind(as.data.frame(st_coordinates(c)),c) [,1:3]
names(c_df) <- c("x_rr", "y_rr", "site_code")
c_df <- c_df %>%
  mutate(
    site = case_when(
      site_code == "CP" ~ "cooks_path",
      site_code == "CB19" ~ "corner_beach",
      site_code == "EP" ~ "easter_point",
      site_code == "GT19" ~ "gnarly_tree",
      site_code == "HS19" ~ "horseshoe",
      site_code == "L119" ~ "lagoon_1",
      site_code == "L219" ~ "lagoon_2",
      site_code == "LH" ~ "lizard_head",
      site_code == "NML" ~ "no_mans_land",
      site_code == "NP" ~ "north_of_paradise",
      site_code == "N1" ~ "north_reef_1",
      site_code == "N2" ~ "north_reef_2",
      site_code == "N319" ~ "north_reef_3",
      site_code == "OS" ~ "osprey",
      site_code == "RS19" ~ "resort",
      site_code == "SI" ~ "south_island",
      site_code == "SE19" ~ "southeast",
      site_code == "TM19" ~ "trimodal",
      site_code == "TB19" ~ "turtle_beach",
      site_code == "VK19" ~ "vickis",
      site_code == "WM" ~ "washing_machine"
    )
  ) %>% 
  filter( site %in% subset_sites)


grid50 <- create_grid(c_df, L = 0.5, "output/grid50_coords_2019_11.csv")
grid1 <- create_grid(c_df, L = 1, "output/grid1_coords_2019_11.csv")
grid2 <- create_grid(c_df, L = 2, "output/grid2_coords_2019_11.csv")
grid4 <- create_grid(c_df, L = 4, "output/grid4_coords_2019_11.csv")


ann <- ann %>% filter(Status != "gone", Species != "soft")
ann_juv <- ann %>% filter(Status == "new")
ann_adu <- ann %>% filter(Status != "new")
counts_2019 <- grid_counts(
  grid_list   = list(grid50, grid1, grid2, grid4),
  square_sizes = c(0.5, 1, 2, 4),
  ann_adu,
  ann_juv,
  output_file = "output/grid_counts_2019.csv"
)

plot_grids_by_site(counts_2019,
                   ann_adu,
                   ann_juv,
                   buffer = 3 ,
                   path_year = 19)


# everything needs to be transformed in meters projections, same as done for maps
ref <- rast("data/dem_reprojected_resampled/TM23_repr.tif")
ann <- st_read("data/lizard_rr_annotations_2023.shp")
ann[is.na(ann$Species), ]$Species <- "Spp"
subset_sites <- c(
  "cooks_path",
  "north_reef_1",
  "gnarly_tree",
  "corner_beach",
  "horseshoe",
  "osprey",
  "resort",
  "south_island",
  "trimodal",
  "turtle_beach",
  "vickis",
  "washing_machine"
)
crs(ann)
ann <- st_transform(ann, crs(ref))
c <- st_read("data/dem_coords_nov2023.shp")
c <- st_transform(c, crs(ref))
c_df <- cbind(as.data.frame(st_coordinates(c)),c) [,1:3]
names(c_df) <- c("x_rr", "y_rr", "site_code")
c_df <- c_df %>%
  mutate(
    site = case_when(
      site_code == "CP" ~ "cooks_path",
      site_code == "CB" ~ "corner_beach",
      site_code == "EP" ~ "easter_point",
      site_code == "GT" ~ "gnarly_tree",
      site_code == "HS" ~ "horseshoe",
      site_code == "L1" ~ "lagoon_1",
      site_code == "L2" ~ "lagoon_2",
      site_code == "LH" ~ "lizard_head",
      site_code == "NML" ~ "no_mans_land",
      site_code == "NP" ~ "north_of_paradise",
      site_code == "N1" ~ "north_reef_1",
      site_code == "N2" ~ "north_reef_2",
      site_code == "N3" ~ "north_reef_3",
      site_code == "OS" ~ "osprey",
      site_code == "RS" ~ "resort",
      site_code == "SI" ~ "south_island",
      site_code == "SE" ~ "southeast",
      site_code == "TM" ~ "trimodal",
      site_code == "TB" ~ "turtle_beach",
      site_code == "VK" ~ "vickis",
      site_code == "WM" ~ "washing_machine"
    )
  ) %>% 
  filter( site %in% subset_sites)


grid50 <- create_grid(c_df, L = 0.5, "output/grid50_coords_2023_11.csv")
grid1 <- create_grid(c_df, L = 1, "output/grid1_coords_2023_11.csv")
grid2 <- create_grid(c_df, L = 2, "output/grid2_coords_2023_11.csv")
grid4 <- create_grid(c_df, L = 4, "output/grid4_coords_2023_11.csv")


ann <- ann %>% filter(Status != "gone", Species != "soft")
ann_juv <- ann %>% filter(ColnyID > 48458)
ann_adu <- ann %>% filter(ColnyID <= 48458)
counts_2023 <- grid_counts(
  grid_list   = list(grid50, grid1, grid2, grid4),
  square_sizes = c(0.5, 1, 2, 4),
  ann_adu,
  ann_juv,
  output_file = "output/grid_counts_2023.csv"
)

plot_grids_by_site(counts_2023,
                   ann_adu,
                   ann_juv,
                   buffer = 3 ,
                   path_year = 23)
