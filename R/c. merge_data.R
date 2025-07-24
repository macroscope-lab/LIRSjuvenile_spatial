### merge counts and complexity 

# upload libraries
library(dplyr)
library(ggplot2)
library(GGally)

# read datasets and merge 2022 #####

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


subset_sites <- c("cooks_path","easter_point","lagoon_2",
                  "north_reef_1","gnarly_tree","corner_beach",
                  "horseshoe", "north_reef_3", "resort", "southeast",
                  "trimodal","turtle_beach")

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

library(gridExtra)
library(grid)
library(cowplot)

pairplot1 <- dt1 %>%
  filter(site_full %in% subset_sites) %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid",
                      lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                                   mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p <- pairplot1

dt1 <- dt1 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt1, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))


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

dt2 <- dt2 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt2, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))


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

dt4 <- dt4 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt4, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))

write.csv(dt1,"data/dt1_2022.csv")
write.csv(dt2,"data/dt2_2022.csv")
write.csv(dt4,"data/dt4_2022.csv")

# read datasets and merge 2019 #####

set <- read.csv("data_raw/rr2019/tiles/settlement_2019-20_LIRS_VB.csv")
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

s <- s %>%
  mutate(site_full = case_when(
    site == "CB19" ~ "corner_beach",
    site == "GT19" ~ "gnarly_tree",
    site == "HS19" ~ "horseshoe",
    site == "L119" ~ "lagoon_1",
    site == "N319" ~ "north_reef_3",
    site == "RS19" ~ "resort",
    site == "SE19" ~ "southeast",
    site == "TM19" ~ "trimodal",
    site == "TB19" ~ "turtle_beach")) %>%
  mutate(site_code = case_when(
    site == "CB19" ~ "CB",
    site == "GT19" ~ "GT",
    site == "HS19" ~ "HS",
    site == "L119" ~ "L1",
    site == "N319" ~ "N3",
    site == "RS19" ~ "RS",
    site == "SE19" ~ "SE",
    site == "TM19" ~ "TM",
    site == "TB19" ~ "TB"))

counts1 <- read.csv("data/grid1_counts19.csv") 
counts2 <- read.csv("data/grid2_counts19.csv") 
counts4 <- read.csv("data/grid4_counts19.csv") 

names(counts1)[1] <- "site_full"
names(counts2)[1] <- "site_full"
names(counts4)[1] <- "site_full"


counts1 <- counts1 %>%
  mutate(site_code = case_when(
    site_full == "cooks_path" ~ "CP",
    site_full == "corner_beach" ~ "CB",
    site_full == "easter_point" ~ "EP",
    site_full == "gnarly_tree" ~ "GT",
    site_full == "horseshoe" ~ "HS",
    site_full == "lagoon_1" ~ "L1",
    site_full == "lagoon_2" ~ "L2",
    site_full == "lizard_head" ~ "LH",
    site_full == "no_mans_land" ~ "NM",
    site_full == "north_of_paradise" ~ "NP",
    site_full == "north_reef_1" ~ "N1",
    site_full == "north_reef_2" ~ "N2",
    site_full == "north_reef_3" ~ "N3",
    site_full == "osprey" ~ "OS",
    site_full == "resort" ~ "RS",
    site_full == "south_island" ~ "SI",
    site_full == "southeast" ~ "SE",
    site_full == "trimodal" ~ "TM",
    site_full == "turtle_beach" ~ "TB",
    site_full == "vickis" ~ "VK",
    site_full == "washing_machine" ~ "WM",
  ))

counts2 <- counts2 %>%
  mutate(site_code = case_when(
    site_full == "cooks_path" ~ "CP",
    site_full == "corner_beach" ~ "CB",
    site_full == "easter_point" ~ "EP",
    site_full == "gnarly_tree" ~ "GT",
    site_full == "horseshoe" ~ "HS",
    site_full == "lagoon_1" ~ "L1",
    site_full == "lagoon_2" ~ "L2",
    site_full == "lizard_head" ~ "LH",
    site_full == "no_mans_land" ~ "NM",
    site_full == "north_of_paradise" ~ "NP",
    site_full == "north_reef_1" ~ "N1",
    site_full == "north_reef_2" ~ "N2",
    site_full == "north_reef_3" ~ "N3",
    site_full == "osprey" ~ "OS",
    site_full == "resort" ~ "RS",
    site_full == "south_island" ~ "SI",
    site_full == "southeast" ~ "SE",
    site_full == "trimodal" ~ "TM",
    site_full == "turtle_beach" ~ "TB",
    site_full == "vickis" ~ "VK",
    site_full == "washing_machine" ~ "WM",
  ))

counts4 <- counts4 %>%
  mutate(site_code = case_when(
    site_full == "cooks_path" ~ "CP",
    site_full == "corner_beach" ~ "CB",
    site_full == "easter_point" ~ "EP",
    site_full == "gnarly_tree" ~ "GT",
    site_full == "horseshoe" ~ "HS",
    site_full == "lagoon_1" ~ "L1",
    site_full == "lagoon_2" ~ "L2",
    site_full == "lizard_head" ~ "LH",
    site_full == "no_mans_land" ~ "NM",
    site_full == "north_of_paradise" ~ "NP",
    site_full == "north_reef_1" ~ "N1",
    site_full == "north_reef_2" ~ "N2",
    site_full == "north_reef_3" ~ "N3",
    site_full == "osprey" ~ "OS",
    site_full == "resort" ~ "RS",
    site_full == "south_island" ~ "SI",
    site_full == "southeast" ~ "SE",
    site_full == "trimodal" ~ "TM",
    site_full == "turtle_beach" ~ "TB",
    site_full == "vickis" ~ "VK",
    site_full == "washing_machine" ~ "WM",
  ))

rdh1 <- read.csv("data/grid1_rdh19.csv") 
rdh2 <- read.csv("data/grid2_rdh19.csv") 
rdh4 <- read.csv("data/grid4_rdh19.csv") 

rdh1 <- rdh1 %>% 
  mutate(site_code = case_when(
  site == "CB19" ~ "CB",
  site == "GT19" ~ "GT",
  site == "HS19" ~ "HS",
  site == "L119" ~ "L1",
  site == "N319" ~ "N3",
  site == "RS19" ~ "RS",
  site == "SE19" ~ "SE",
  site == "TM19" ~ "TM",
  site == "TB19" ~ "TB"))

rdh2 <- rdh2 %>% 
  mutate(site_code = case_when(
    site == "CB19" ~ "CB",
    site == "GT19" ~ "GT",
    site == "HS19" ~ "HS",
    site == "L119" ~ "L1",
    site == "N319" ~ "N3",
    site == "RS19" ~ "RS",
    site == "SE19" ~ "SE",
    site == "TM19" ~ "TM",
    site == "TB19" ~ "TB"))

rdh4 <- rdh4 %>% 
  mutate(site_code = case_when(
    site == "CB19" ~ "CB",
    site == "GT19" ~ "GT",
    site == "HS19" ~ "HS",
    site == "L119" ~ "L1",
    site == "N319" ~ "N3",
    site == "RS19" ~ "RS",
    site == "SE19" ~ "SE",
    site == "TM19" ~ "TM",
    site == "TB19" ~ "TB"))

dt1 <- merge(counts1,rdh1, by=c("site_code","rep"))
dt2 <- merge(counts2,rdh2, by=c("site_code","rep"))
dt4 <- merge(counts4,rdh4, by=c("site_code","rep"))

dt1 <- dt1[,!(names(dt1) %in% c("X", "site_full.y", "site_code.y")) ]
dt2 <- dt2[,!(names(dt2) %in% c("X", "site_full.y", "site_code.y")) ]
dt4 <- dt4[,!(names(dt4) %in% c("X",  "site_full.y", "site_code.y")) ]

dt1 <- merge(dt1,s, by="site")
dt2 <- merge(dt2,s, by="site")
dt4 <- merge(dt4,s, by="site")

dt1 <- dt1[,!(names(dt1) %in% c("X", "site_full.y", "site_code.y")) ]
dt2 <- dt2[,!(names(dt2) %in% c("X", "site_full.y", "site_code.y")) ]
dt4 <- dt4[,!(names(dt4) %in% c("X",  "site_full.y", "site_code.y")) ]

names(dt1)[2] <- "site_code"
names(dt2)[2] <- "site_code"
names(dt4)[2] <- "site_code"

names(dt1)[4] <- "site_full"
names(dt2)[4] <- "site_full"
names(dt4)[4] <- "site_full"

subset_sites19 <- c("lagoon_1",
                  "north_reef_3","gnarly_tree","corner_beach",
                  "horseshoe", "resort", "southeast",
                  "trimodal","turtle_beach")

dt1 <- dt1 %>%
  filter(site_full %in% subset_sites19) %>%
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

library(gridExtra)
library(grid)
library(cowplot)

pairplot1_19 <- dt1[dt1$year == 2019,] %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid - 19",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p19 <- pairplot1_19


pairplot1_22 <- dt1[dt1$year == 2022,] %>%
  filter(site != "L2") %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid - 2022",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p22 <- pairplot1_22


dt1 <- dt1 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt1, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))


ggsave(file="figures/Pairplot1_19.png", plot=pairplot1, width=7, height=6)

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

dt2 <- dt2 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt2, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))


ggsave(file="figures/Pairplot2_19.png", plot=pairplot2, width=7, height=6)




pairplot2_19 <- dt2[dt2$year == 2019,] %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 2m grid - 19",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p219 <- pairplot2_19


pairplot2_22 <- dt2[dt2$year == 2022,] %>%
  filter(site != "L2") %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 2m grid - 2022",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p222 <- pairplot2_22


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

ggsave(file="figures/Pairplot4_19.png", plot=pairplot4, width=7, height=6)

dt4 <- dt4 %>%
  mutate(site = reorder(site, R, FUN = median)) %>%
  group_by(site) %>%
  mutate (medianD = median(D),
          medianR = median(R))

dummy_plot <-
  ggplot(dt4, aes(x = R, y = D, color = site)) + 
  geom_point() + scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme(legend.position = "right")+
  theme_minimal()

legend <- cowplot::get_legend(dummy_plot)

pairplot4_19 <- dt4_19 %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 4m grid - 19",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p419 <- pairplot4_19


pairplot4_22 <- dt4_22 %>%
  filter(site != "L2") %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 4m grid - 2022",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p422 <- pairplot4_22


# Display with legend using cowplot
cowplot::plot_grid(p, legend, rel_widths = c(4, 0.4))

dt1_19 <- dt1
dt2_19 <- dt2
dt4_19 <- dt4

write.csv(dt1_19,"data/dt1_2019.csv")
write.csv(dt2_19,"data/dt2_2019.csv")
write.csv(dt4_19,"data/dt4_2019.csv")

### merge years ####
dt1_22 <- read.csv("data/dt1_2022.csv")
dt2_22 <- read.csv("data/dt2_2022.csv")
dt4_22 <- read.csv("data/dt4_2022.csv")

dt1_19 <- read.csv("data/dt1_2019.csv")
dt2_19 <- read.csv("data/dt2_2019.csv")
dt4_19 <- read.csv("data/dt4_2019.csv")

dt1_19$year <- 2019
dt2_19$year <- 2019
dt4_19$year <- 2019

dt1_22$year <- 2022
dt2_22$year <- 2022
dt4_22$year <- 2022

names(dt1_22)
names(dt1_19)

dt1_22 <- dt1_22[,c(2:24)]
dt1_19 <- dt1_19[,c(3:25)]
names(dt1_19)[1] <- "site"
column_order <- names(dt1_19)
dt1_22 <- dt1_22[, c(column_order)]
names(dt1_19)
names(dt1_22)

dt1 <- rbind (dt1_19, dt1_22)
write.csv(dt1, "dt1.csv")

names(dt2_22)
names(dt2_19)

dt2_22 <- dt2_22[,c(2:20)]
dt2_19 <- dt2_19[,c(3:21)]
names(dt2_19)[1] <- "site"
column_order <- names(dt2_19)
dt2_22 <- dt2_22[, c(column_order)]
names(dt2_19)
names(dt2_22)

dt2 <- rbind (dt2_19, dt2_22)
write.csv(dt2, "dt2.csv")

names(dt4_22)
names(dt4_19)

dt4_22 <- dt4_22[,c(2:13,15:20)]
dt4_19 <- dt4_19[,c(3:20)]
names(dt4_19)[1] <- "site"
column_order <- names(dt4_19)
dt4_22 <- dt4_22[, c(column_order)]
names(dt4_19)
names(dt4_22)

dt4 <- rbind (dt4_19, dt4_22)
write.csv(dt4, "dt4.csv")


pairplot1 <- dt1 %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site","year") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 1m grid",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p1 <- pairplot1
ggsave(file="figures/Pairplot1_all.png", plot=pairplot1, width=9, height=8)


pairplot2 <- dt2 %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site","year") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 2m grid",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p2 <- pairplot2
ggsave(file="figures/Pairplot2_all.png", plot=pairplot2, width=9, height=8)

dt4_19
pairplot4_19 <- dt4_19 %>%
  select("n_adu", "n_juv","mean_sett","R","D","H","site","year") %>%
  mutate(site = reorder(site, R, FUN = median))%>%
  ggpairs(title = "Pair plot - 4m grid",
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=site)), 
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
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))
p419 <- pairplot4_19
ggsave(file="figures/Pairplot4_all.png", plot=pairplot4, width=9, height=8)



library(ggrain)
dt1 %>% filter(site != "L2") %>%
ggplot(aes(1, n_adu, fill = site)) +
  geom_rain(alpha = .3, rain.side = 'l',
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .1, width = 0.1), width = 0.1
            )) +
  theme_bw() +
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_grid(year~site)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


ggplot(dt2,aes(1, n_adu, fill = site)) +
  geom_rain(alpha = .3, rain.side = 'l',
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .1, width = 0.1), width = 0.1
            )) +
  theme_bw() +
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_grid(year~site)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())




dt1 %>% filter(site != "L2") %>%
  ggplot(aes(1, n_juv, fill = site)) +
  geom_rain(alpha = .3, rain.side = 'l',
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .1, width = 0.1), width = 0.1
            )) +
  theme_bw() +
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_grid(year~site)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())





dt1 %>% select(n_adu, year, site) %>% 
  ggplot() +
  geom_boxplot(aes(n_adu))
