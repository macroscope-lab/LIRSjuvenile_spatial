### merge counts and complexity 

# upload libraries
library(dplyr)
library(ggplot2)

# read datasets and merge
counts1 <- read.csv("data/grid1_counts.csv") 
counts2 <- read.csv("data/grid2_counts.csv") 
counts4 <- read.csv("data/grid4_counts.csv") 

counts1$site_full <- counts1$site
counts2$site_full <- counts2$site
counts4$site_full <- counts4$site

counts1 <- counts1 %>%
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

counts2 <- counts2 %>%
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

counts4 <- counts4 %>%
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

rdh1 <- read.csv("data/grid1_rdh.csv") 
rdh2 <- read.csv("data/grid2_rdh.csv") 
rdh4 <- read.csv("data/grid4_rdh.csv") 

dt1 <- merge(counts1,rdh1, by=c("site","rep"))
dt2 <- merge(counts2,rdh2, by=c("site","rep"))
dt4 <- merge(counts4,rdh4, by=c("site","rep"))

dt1 <- dt1[,!(names(dt1) %in% c("X", "X.1")) ]
dt2 <- dt2[,!(names(dt2) %in% c("X", "X.1")) ]
dt4 <- dt4[,!(names(dt4) %in% c("X", "X.1")) ]

# plot variables






