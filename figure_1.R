# Figure 1 : Map and broad-scale population structure

library(ggpubr)
library(sf)
library(ozmaps)
library(grid)
source("constants.R")
library(ggrepel)
library(tidyverse)

load("cache/sites_data.rdata")
run12_meta <- read_rds("cache/run12_meta.rds")

adjacent_dispersal <- read_sf("dispersal/ALL_TOGETHER_footprint.kml",layer = "Line Features")
maggie_dispersal <- read_sf("dispersal/ALL_TOGETHER_footprint.kml",layer = "Area Features")


gbr_coast <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "GBR_coastlineONLY")
#gbr_tc <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "GBR_towncities")
gbr_features <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")

#sf_use_s2(FALSE) # the input geometries are bad, and wrong spherically
x1 <- st_make_valid(gbr_features)

xmn=146.0
ymn=-19.3
xmx=146.9
ymx=-18.35



# xmn=146.2
# ymn=-19.4
# xmx=147.6
# ymx=-17.5

gbr_ft_cropped <- st_crop(x1, xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx) #%>% filter(LEVEL_1=="Island")
gbr_coast_cropped <- st_crop(gbr_coast,xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx)

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sites_sf <- st_as_sf(sites_data,coords = c("lon", "lat"),crs = projcrs) %>% 
  mutate(name_label = paste(name," (",pop,")",sep=""))

pmap <-   ggplot() + 
#  geom_sf(data = adjacent_dispersal, color = "blue", fill="transparent", size=0.5) +
#  geom_sf(data = maggie_dispersal, color = "orange", fill="transparent",size=0.5) +    
  geom_sf(data = gbr_ft_cropped %>% filter(FEAT_NAME %in% c("Reef","Island")), aes(fill=FEAT_NAME)) +
  geom_sf(data = gbr_coast_cropped) +
  geom_sf(data = sites_sf,aes(color=pop),size=2) + 
  scale_fill_manual(values = c("Reef"="945200","Island"="grey90")) +
  scale_color_manual(values = location_colors) +
  ggrepel::geom_text_repel(data = sites_sf, mapping = aes(label=name_label, geometry=geometry), stat = "sf_coordinates",size=2, force = 10) + 
  theme_bw() +
  xlab("") + ylab("")

padmix <- run12_meta %>% 
  mutate(pop_order = location_order[pop]) %>% 
  ggplot(aes(y = ID, x = p)) +
  scale_fill_manual(values = cluster_colors, labels = cluster_names) +
  geom_col(aes(fill=cluster),position="stack", stat = "identity",width = 1,linewidth=0) +
  facet_grid(rows=vars(reorder(pop,pop_order)),scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom", axis.text.x = element_text(size=6)) + 
  theme(strip.text = element_text(size=6)) +
  ylab("") + xlab("Admixture Proportion") + labs(fill=NULL)


g <- ggplot_gtable(ggplot_build(padmix + theme(legend.position="none")))

strip <- which(grepl('strip-', g$layout$name))

fills <- location_colors
k <- 1
for (i in strip) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)


library(cowplot)

plot_grid(
  pmap + theme(legend.position="none"),
  g,
  rel_widths = c(3, 0.5))

ggsave(filename = "figures/Figure1.png",width = 174,height = 110, units = "mm")
ggsave(filename = "figures/Figure1.svg",width = 174,height = 110,units = "mm")


plot_grid(
  pmap + theme(legend.position="none"),
  padmix + theme(legend.position = "none"),
  rel_widths = c(3, 1))
