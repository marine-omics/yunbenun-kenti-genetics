# Figure S1: Hydrodynamic modelling


library(ggpubr)
library(sf)
library(ozmaps)
source("constants.R")
library(ggrepel)
library(tidyverse)

load("cache/sites_data.rdata")
sea_dist <- read_rds("cache/sea_dist.rds")

c("BRR","EP","GB","HaR","HB","HFB","JB","KR","LPB","MB","MR","PB","SO","WB","WP")

maggie_kmls <- c("HFB"="Huntingfield_Bay.kml","WB"="Wilson_Bay.kml","MB"="Maud_Bay.kml","HB"="Horseshoe_Bay.kml","GB"="Geoffrey_Bay.kml","PB"="Picnic_Bay.kml","MR"="Middle_Reef.kml")
palms_kmls <- c("BRR"="Bramble_Reef.kml","LPB"="Orpheus_Island.kml","HaR"="Havannah_Reef.kml","JB"="John_Brewer_Reef.kml","KR"="Keeper_Reef.kml","WP"="Pelorus_Reef.kml")


kml_files <- paste("dispersal",c(maggie_kmls,palms_kmls),sep="/")
names(kml_files) <- c(names(maggie_kmls),names(palms_kmls))

disp_data <- map(kml_files,st_read) %>% bind_rows(.id = "pop")

disp_sf_sites <- disp_data %>% 
  left_join(sites_data,by="pop") %>% 
  mutate()

disp_adj <- disp_sf_sites %>% 
  filter(pop_group=="ADJ") %>% 
  st_combine()

disp_mi <- disp_sf_sites %>% 
  filter(pop_group=="MI") %>% 
  st_union()

gbr_coast <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "GBR_coastlineONLY")
gbr_features <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")

sf_use_s2(FALSE) # the input geometries are bad, and wrong spherically
x1 <- st_make_valid(gbr_features)

xmn=146.2
ymn=-19.5
xmx=147.5
ymx=-17.52

gbr_ft_cropped <- st_crop(x1, xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx) #%>% filter(LEVEL_1=="Island")
gbr_coast_cropped <- st_crop(gbr_coast,xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx)

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sites_sf <- st_as_sf(sites_data,coords = c("lon", "lat"),crs = projcrs) %>% 
  mutate(name_label = paste(name," (",pop,")",sep=""))

ggplot() + 
  geom_sf(data = gbr_ft_cropped %>% filter(FEAT_NAME %in% c("Reef","Island")), aes(fill=FEAT_NAME)) +
  geom_sf(data = gbr_coast_cropped) +
  geom_sf(data = disp_sf_sites ,aes(color=pop),fill="transparent", size = 1) +
  geom_sf(data = sites_sf ,aes(color=pop),size=2) + 
  ggrepel::geom_text_repel(data = sites_sf %>% filter(pop_group=="ADJ") , mapping = aes(label=name_label, geometry=geometry), stat = "sf_coordinates",size=2, max.overlaps = 50) + 
  scale_color_manual(values = location_colors) +
  theme_bw() +
  guides(color="none") + 
  labs(fill=NULL) +
  theme(legend.position = "bottom") +
  xlab("") + ylab("")

ggsave(filename = "figures/FigureS1.png",width = 10,height = 8)

