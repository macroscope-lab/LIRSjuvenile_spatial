### merge counts and complexity 

# upload libraries
library(dplyr)
library(ggplot2)
library(GGally)

# read datasets and merge 2022 #####

set19 <- read.csv("data_raw/rr2019/tiles/settlement_2019-20_LIRS_VB.csv")
s19 <- set19 %>% 
  group_by(site_abbreviation, tile_number) %>%
  mutate(tile_total = sum(total)) %>%
  dplyr::select(site_abbreviation, tile_number, tile_total)%>%
  unique() %>%
  ungroup() %>%
  group_by(site_abbreviation) %>%
  mutate(mean_settlers = mean(tile_total)) %>%
  dplyr::select(site_abbreviation, mean_settlers) %>%
  unique()
names(s19) <- c("site","mean_sett") 
s19$site <- substr(s19$site, 1,2)
counts19 <- read.csv("output/grid_counts_2019.csv") 
names(counts19)
counts19 <- counts19 %>%
  mutate(site = case_when(
    site == "cooks_path" ~ "CP",
    site == "corner_beach" ~ "CB",
    site == "easter_point" ~ "EP",
    site == "gnarly_tree" ~ "GT",
    site == "horseshoe" ~ "HS",
    site == "lagoon_1" ~ "L1",
    site == "lagoon_2" ~ "L2",
    site == "lizard_head" ~ "LH",
    site == "no_mans_land" ~ "NM",
    site == "north_of_paradise" ~ "NP",
    site == "north_reef_1" ~ "N1",
    site == "north_reef_2" ~ "N2",
    site == "north_reef_3" ~ "N3",
    site == "osprey" ~ "OS",
    site == "resort" ~ "RS",
    site == "south_island" ~ "SI",
    site == "southeast" ~ "SE",
    site == "trimodal" ~ "TM",
    site == "turtle_beach" ~ "TB",
    site == "vickis" ~ "VK",
    site == "washing_machine" ~ "WM",
  ))
geom19_50 <- read.csv("data/grid50_rdh19.csv")
geom19_50$scale <- 0.5
geom19_1 <- read.csv("data/grid1_rdh19.csv")
geom19_1$scale <- 1
geom19_2 <- read.csv("data/grid2_rdh19.csv")
geom19_2$scale <- 2
geom19_4 <- read.csv("data/grid4_rdh19.csv")
geom19_4$scale <- 4

geom19 <- rbind(geom19_50, geom19_1,
                geom19_2, geom19_4)
geom19$site <- substr(geom19$site, 1,2)

merged19 <- geom19 %>%
  full_join(counts19, by = c("site", "rep", "scale"))

merged19 <- merged19 %>% 
  left_join(s19, by = "site") %>%
  mutate (year = 2019)

# 2023
set23 <- read.csv("data_raw/rr2023/tiles/Settlement_LIRS_2022-2023_with_bluelight.csv")
s23 <- set23 %>% 
  group_by(site_abbreviation, tile_number) %>%
  mutate(tile_total = sum(total)) %>%
  dplyr::select(site_abbreviation, tile_number, tile_total)%>%
  unique() %>%
  ungroup() %>%
  group_by(site_abbreviation) %>%
  mutate(mean_settlers = mean(tile_total)) %>%
  dplyr::select(site_abbreviation, mean_settlers) %>%
  unique()
names(s23) <- c("site","mean_sett") 

counts23 <- read.csv("output/grid_counts_2023.csv") 
names(counts23)
counts23 <- counts23 %>% filter(site != "corner_beach") %>%
  mutate(site = case_when(
    site == "cooks_path" ~ "CP",
    site == "corner_beach" ~ "CB",
    site == "easter_point" ~ "EP",
    site == "gnarly_tree" ~ "GT",
    site == "horseshoe" ~ "HS",
    site == "lagoon_1" ~ "L1",
    site == "lagoon_2" ~ "L2",
    site == "lizard_head" ~ "LH",
    site == "no_mans_land" ~ "NM",
    site == "north_of_paradise" ~ "NP",
    site == "north_reef_1" ~ "N1",
    site == "north_reef_2" ~ "N2",
    site == "north_reef_3" ~ "N3",
    site == "osprey" ~ "OS",
    site == "resort" ~ "RS",
    site == "south_island" ~ "SI",
    site == "southeast" ~ "SE",
    site == "trimodal" ~ "TM",
    site == "turtle_beach" ~ "TB",
    site == "vickis" ~ "VK",
    site == "washing_machine" ~ "WM",
  ))
  

geom23_50 <- read.csv("data/grid50_rdh23.csv")
geom23_50N3 <- read.csv("data/grid50_rdh23N3.csv")
geom23_50 <- rbind(geom23_50,geom23_50N3, geom23_50CB)
geom23_50$scale <- 0.5
geom23_1 <- read.csv("data/grid1_rdh23.csv")
geom23_1N3 <- read.csv("data/grid1_2023_rdhN3.csv")
geom23_1 <- rbind(geom23_1,geom23_1N3, geom23_1CB)
geom23_1$scale <- 1
geom23_2 <- read.csv("data/grid2_rdh23.csv")
geom23_2N3 <- read.csv("data/grid2_2023_rdhN3.csv")
geom23_2 <- rbind(geom23_2,geom23_2N3, geom23_2CB)
geom23_2$scale <- 2
geom23_4 <- read.csv("data/grid4_rdh23.csv")
geom23_4N3 <- read.csv("data/grid4_2023_rdhN3.csv")
geom23_4 <- rbind(geom23_4,geom23_4N3, geom23_4CB)
geom23_4$scale <- 4

geom23 <- rbind(geom23_50, geom23_1,
                geom23_2, geom23_4)

merged23 <- geom23 %>%
  full_join(counts23, by = c("site", "rep", "scale"))

merged23 <- merged23 %>% 
  left_join(s23, by = "site") %>%
  mutate (year = 2023)

dt <- rbind(merged19, merged23)
dt$year <- as.factor(dt$year)
dt$site <- as.factor(dt$site)
dt <- dt %>%
  mutate(R_log10 = log10(R),
         H_log10 = log10(H),
         n_adu_log = log(n_adu+1),
         n_juv_log = log(n_juv+1),
         mean_sett_log = log(mean_sett))

dt50 <- dt %>% filter(scale == 0.5) %>%
  mutate(Ds = as.numeric(scale(D)),
         n_adu_s = scale(n_adu))
write.csv(dt50,"output/dt50.csv")

dt1 <- dt %>% filter(scale == 1) %>%
  mutate(Ds = as.numeric(scale(D)),
         n_adu_s = scale(n_adu))
write.csv(dt1, "output/dt1.csv")

dt2 <- dt %>% filter(scale == 2) %>%
  mutate(Ds = as.numeric(scale(D)),
         n_adu_s = scale(n_adu))
write.csv(dt2, "output/dt2.csv")

dt4 <- dt %>% filter(scale ==4) %>%
  mutate(Ds = as.numeric(scale(D)),
         n_adu_s = scale(n_adu))
write.csv(dt4, "output/dt4.csv")

#### plot variables

library(gridExtra)
library(grid)
library(cowplot)

pp50 <- dt50 %>%
  dplyr::select("n_adu_log", "n_juv_log","R_log10","Ds","H_log10","mean_sett_log", "site","year") %>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(title = "Pair plot - 50cm grid",
                      lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size=0.5),  
                                   mapping = ggplot2::aes(colour=site)), 
                      diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                                  continuous = GGally::wrap("densityDiag", alpha=0.5)), 
                      upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                                   continuous = GGally::wrap("cor", size=4, alignPercent=0.8)),
                      labeller = "label_parsed", 
                     cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pp50

pdf("figures/pairplot50.pdf", height = 7, width = 12)
print(pp50)
dev.off()


pp1 <- dt1 %>%
  dplyr::select("n_adu_log", "n_juv_log","R_log10","Ds","H_log10","mean_sett_log", "site","year") %>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid",
          lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
          diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                      continuous = GGally::wrap("densityDiag", alpha=0.5)), 
          upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8)),
          labeller = "label_parsed", 
          cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pp1

pdf("figures/pairplot1.pdf", height = 7, width = 12)
print(pp1)
dev.off()


pp2 <- dt2 %>%
  dplyr::select("n_adu_log", "n_juv_log","R_log10","Ds","H_log10","mean_sett_log", "site","year") %>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(title = "Pair plot - 2m grid",
          lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
          diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                      continuous = GGally::wrap("densityDiag", alpha=0.5)), 
          upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8)),
          labeller = "label_parsed", 
          cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pp2

pdf("figures/pairplot2.pdf", height = 7, width = 12)
print(pp2)
dev.off()

pp4 <- dt4 %>%
  dplyr::select("n_adu_log", "n_juv_log","R_log10","Ds","H_log10","mean_sett_log", "site","year") %>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(title = "Pair plot - 4m grid",
          lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
          diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                      continuous = GGally::wrap("densityDiag", alpha=0.5)), 
          upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8)),
          labeller = "label_parsed", 
          cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pp4


pdf("figures/pairplot4.pdf", height = 7, width = 12)
print(pp4)
dev.off()

