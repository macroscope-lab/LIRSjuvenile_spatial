### merge counts and complexity 

# upload libraries
library(dplyr)
library(ggplot2)
library(GGally)

# read datasets and merge

set <- read.csv("data_raw/Settlement_LIRS_2021-2022_VB_MJB_JC_MW_DW.csv")
s <- set %>% 
  group_by(site_abbreviation, tile_number) %>%
  mutate(tile_total = sum(total)) %>%
  select(site_abbreviation, tile_number, tile_total)%>%
  unique() %>%
  ungroup() %>%
  group_by(site_abbreviation) %>%
  mutate(mean_settlers = mean(tile_total)) %>%
  select(site_abbreviation, mean_settlers) %>%
  unique()
s <- s[1:20,]
names(s) <- c("site","mean_sett") 

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

dt1 <- merge(dt1,s, by="site")
dt2 <- merge(dt2,s, by="site")
dt4 <- merge(dt4,s, by="site")
  
dt1 <- dt1 %>%
  filter(site_full %in% subset_sites) %>%
  mutate(
  n_adu_log = log(n_adu+1),
  mean_sett_s = as.numeric(scale(mean_sett)),
  R_log10 = log10(R),
  D_s = as.numeric(scale(D))) %>%
  group_by(site)%>%
  mutate(xs = mean(xc),
         ys = mean(yc)) %>%
  ungroup() 


dt2 <-  dt2 %>%
  filter(site_full %in% subset_sites) %>%
  group_by(site)%>%
  mutate(xs = mean(xc),
         ys = mean(yc)) %>%
  ungroup()


dt4 <- dt4 %>%
  filter(site_full %in% subset_sites) %>%
  group_by(site)%>%
  mutate(xs = mean(xc),
         ys = mean(yc)) %>%
  ungroup()


#### plot variables

subset_sites <- c("cooks_path","easter_point","lagoon_2",
                  "north_reef_1","gnarly_tree","corner_beach",
                  "horseshoe", "north_reef_3", "resort", "southeast",
                  "trimodal","turtle_beach")

pairplot1 <- dt1 %>%
  filter(site_full %in% subset_sites) %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid",
                      lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=site)), 
                      diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                                  continuous = wrap("densityDiag", alpha=0.5)), 
                      upper = list(combo = wrap("box_no_facet", alpha=0.5),
                                   continuous = wrap("cor", size=4, alignPercent=0.8)),
                      #labeller = "label_parsed", 
                     cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pairplot1

ggsave(file="figures/Pairplot1.png", plot=pairplot1, width=7, height=6)

pairplot2 <- dt2 %>%
  filter(site_full %in% subset_sites) %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 2m grid",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=site)), 
          diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                      continuous = wrap("densityDiag", alpha=0.5)), 
          upper = list(combo = wrap("box_no_facet", alpha=0.5),
                       continuous = wrap("cor", size=4, alignPercent=0.8)),
          #labeller = "label_parsed", 
          cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pairplot2

ggsave(file="figures/Pairplot2.png", plot=pairplot2, width=7, height=6)

pairplot4 <- dt4 %>%
  filter(site_full %in% subset_sites) %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 4m grid",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  mapping = ggplot2::aes(colour=site)), 
          diag = list(discrete="barDiag",  mapping = ggplot2::aes(colour=site),
                      continuous = wrap("densityDiag", alpha=0.5)), 
          upper = list(combo = wrap("box_no_facet", alpha=0.5),
                       continuous = wrap("cor", size=4, alignPercent=0.8)),
          #labeller = "label_parsed", 
          cardinality_threshold = 18) + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
pairplot4

ggsave(file="figures/Pairplot4.png", plot=pairplot4, width=7, height=6)



