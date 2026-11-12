library(ggplot2)
library(GGally)
library(dplyr)
library(tidyverse)


## pairplot
paths <- c("output/dt1.csv", "output/dt50.csv", "output/dt2.csv")

dt <- read.csv(paths[1])
dt$site <- as.factor(dt$site)
dt$year <- as.factor(dt$year)
dt$n_adu_log_c     <- scale(dt$n_adu_log, scale = FALSE)
dt$mean_sett_log_c <- scale(dt$mean_sett_log, scale = FALSE)
dt$R_log10_c       <- scale(dt$R_log10, scale = FALSE)
dt$Dc            <- scale(dt$D, scale = FALSE)
dt$n_adu_log_s     <- scale(dt$n_adu_log)
dt$mean_sett_log_s <- scale(dt$mean_sett_log)
dt$R_log10_s       <- scale(dt$R_log10)
dt$Ds            <- scale(dt$D)

dt %>%
  dplyr::select("n_juv_log","n_adu_log", "R_log10","D","mean_sett_log", "site","year") %>%
  mutate(sampling_event = paste0(site, year),
         year = as.factor(year))%>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(columns = 1:5,
    lower = list(continuous = GGally::wrap("smooth", alpha = 0.3, size=0.5),  
                       mapping = ggplot2::aes(colour=sampling_event)
                       ), 
          diag = list(discrete="barDiag",  
                      mapping = ggplot2::aes(colour=sampling_event),
                      continuous = GGally::wrap("densityDiag", alpha=0.5)
                      ), 
          upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8),  
                       mapping = ggplot2::aes(fill=sampling_event)),
          labeller = as_labeller(c(
            n_juv_log     = "juvenile \nabundance \n(log)",
            n_adu_log     = "adult \nabundance\n(log)",
            R_log10       = "surface \nrugosity\n(log₁₀)",
            D             = "fractal \ndimension",
            mean_sett_log = "average \nsettlement\n(log)",
            site          = "site",
            year          = "year",
            sampling_event = "sampling \nevent"
          )),
          cardinality_threshold = 21) + 
  
  scale_colour_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  scale_fill_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 10,
      face = "plain"
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 8
    ),
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.title = element_text(
      colour = "black",
      size = 10
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    )
  )

ggsave("figures/pp1m_1.png",
  width = 17,
  height = 17,
  units = "cm",
  dpi = 300,
  limitsize = TRUE)

dt %>%
  dplyr::select("n_juv_log","n_adu_log", "R_log10","D","mean_sett_log", "site","year") %>%
  mutate(sampling_event = paste0(site, year),
         year = as.factor(year))%>%
  mutate(site = reorder(site, R_log10, FUN = median))%>%
  ggpairs(columns = 1:8,
          lower = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8),  
                       mapping = ggplot2::aes(fill=sampling_event)),
          diag = list(discrete="barDiag",  
                      mapping = ggplot2::aes(colour=sampling_event),
                      continuous = GGally::wrap("densityDiag", alpha=0.5)
          ), 
          upper = list(combo = GGally::wrap("box_no_facet", alpha=0.5),
                       continuous = GGally::wrap("cor", size=4, alignPercent=0.8),  
                       mapping = ggplot2::aes(fill=sampling_event)),
          labeller = as_labeller(c(
            n_juv_log     = "juvenile \nabundance \n(log)",
            n_adu_log     = "adult \nabundance\n(log)",
            R_log10       = "surface \nrugosity\n(log₁₀)",
            D             = "fractal \ndimension",
            mean_sett_log = "Average \nsettlement\n(log)",
            site          = "site",
            year          = "year",
            sampling_event = "sampling \nevent"
          )),
          cardinality_threshold = 21) + 
  
  scale_colour_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  scale_fill_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 10,
      face = "plain"
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 7
    ),
    axis.text.x = element_text(
      angle = 270, vjust = 0.5, hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.title = element_text(
      colour = "black",
      size = 10
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    )
  )
ggsave("figures/pp1m_2.png",
       width = 17,
       height = 35,
       units = "cm",
       dpi = 300,
       limitsize = TRUE)
  
  
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = NULL, colour = "black"),
    strip.text = element_text(size = 12))


library(dplyr)
library(GGally)
library(ggplot2)

# Variables to show in the pair plot
plot_vars <- c(
  "n_juv",
  "n_adu",
  "R",
  "D",
  "mean_sett",
  "site",
  "year"
)

dt_plot <- dt %>%
  dplyr::select(all_of(plot_vars)) %>%
  mutate(
    sampling_event = paste0(site, year),
    year = factor(year),
    site = reorder(site, R, FUN = median)
  )

dt_plot %>%
  ggpairs(
    columns = 1:7,
    title = NULL,
    
    # Continuous-continuous panels
    lower = list(
      continuous = GGally::wrap(
        "smooth",
        alpha = 0.25,
        size = 0.5,
        mapping = aes(colour = sampling_event)
      )
    ),
    
    # Diagonal
    diag = list(
      discrete = GGally::wrap(
        "barDiag",
        mapping = aes(fill = sampling_event)
      ),
      continuous = GGally::wrap(
        "densityDiag",
        alpha = 0.5,
        mapping = aes(fill = sampling_event)
      )
    ),
    
    # Upper panels
    upper = list(
      combo = GGally::wrap(
        "box_no_facet",
        alpha = 0.55,
        mapping = aes(fill = sampling_event)
      ),
      continuous = GGally::wrap(
        "cor",
        size = 4,
        alignPercent = 0.8
      )
    ),
    
    # Rename variables in the grid
    labeller = as_labeller(c(
      n_juv_log     = "juvenile \nabundance \n(log)",
      n_adu_log     = "adult \nabundance\n(log)",
      R_log10       = "surface \nrugosity\n(log₁₀)",
      D             = "fractal \ndimension",
      mean_sett_log = "average \nsettlement\n(log)",
      site          = "site",
      year          = "year",
      sampling_event = "sampling \nevent"
    )),
    
    cardinality_threshold = 21
  ) +
  
  scale_colour_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  scale_fill_viridis_d(
    name = "Sampling event"#,
    #option = "D"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      face = "plain"
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      colour = "black",
      size = 10
    ),
    axis.title = element_text(
      colour = "black",
      size = 10
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    )
  )


library (ggridges)

vars <- c("n_juv_log", "n_adu_log", "R_log10", "D")

dt_long <- dt %>%
select(site, year, all_of(vars)) %>%
pivot_longer(
cols = all_of(vars),
names_to = "variable",
values_to = "value"
) %>%
filter(!is.na(value)) %>%
mutate(
site = factor(site),
year = factor(year, levels = c("2019", "2023")),
variable = factor(
variable,
levels = vars,
labels = c(
"Juvenile\nabundance\n(log)",
"Adult\nabundance\n(log)",
"Surface\nrugosity\n(log₁₀)",
"Fractal\ndimension"
)
)
)

ggplot(dt_long,aes(x = value,y = site,fill = year,colour = year)) +
geom_density_ridges(aes(height = after_stat(density)),stat = "density",position = "identity",scale = 1.15,
alpha = 0.30,linewidth = 0.5,rel_min_height = 0) +
# Observation ticks beneath each density
geom_segment(aes(x = value,xend = value,y = as.numeric(site),yend = as.numeric(site) - 0.15),
linewidth = 0.25,alpha = 0.35) +
scale_fill_manual(values = c("2019" = "#BF812D","2023" = "#35978F")) +
scale_colour_manual(values = c("2019" = "#8C510A","2023" = "#01665E")) +
facet_wrap(~ variable,scales = "free_x",nrow = 1,)+
labs(x = NULL,y = "Site",fill = "Year",colour = "Year") +
theme_classic(base_size = 10) +
theme(legend.position = "bottom",strip.background = element_blank(),
strip.text = element_text(face = "bold",size = 10),
axis.text.y = element_text(),
panel.spacing = unit(1.1, "lines"))


ggsave("figures/fig2_rawdata_2m.png",
width = 14,
height = 14,
units = "cm",
dpi = 300,
limitsize = TRUE)




ggplot(
  dt,
  aes(
    x = value,
    y = site,
    fill = year,
    colour = year
  )
) +
  geom_density_ridges(
    aes(height = after_stat(density)),
    stat = "density",
    position = "identity",
    scale = 1.15,
    alpha = 0.30,
    linewidth = 0.5,
    rel_min_height = 0
  ) +
  # Observation ticks beneath each density
  geom_segment(
    aes(
      x = value,
      xend = value,
      y = as.numeric(site),
      yend = as.numeric(site) - 0.15
    ),
    linewidth = 0.25,
    alpha = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "2019" = "#BF812D",
      "2023" = "#35978F"
    )
  ) +
  scale_colour_manual(
    values = c(
      "2019" = "#8C510A",
      "2023" = "#01665E"
    )
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_x",
    nrow = 1,
  )+
  labs(
    x = NULL,
    y = "Site",
    fill = "Year",
    colour = "Year"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(),
    panel.spacing = unit(1.1, "lines"))


ggsave("figures/fig2_rawdata.png",
       width = 15,
       height = 12,
       units = "cm",
       dpi = 300,
       limitsize = TRUE)




library(ggplot2)
library(GGally)
library(dplyr)
paths <- c("output/dt1.csv", "output/dt50.csv", "output/dt2.csv")
dt <- read.csv(paths[1])
dt$site <- as.factor(dt$site)
dt$year <- as.factor(dt$year)
dt$n_adu_log_c     <- scale(dt$n_adu_log, scale = FALSE)
dt$mean_sett_log_c <- scale(dt$mean_sett_log, scale = FALSE)

