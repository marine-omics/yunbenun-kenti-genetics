# Figure 2: Fine scale population structure


library(ggpubr)
library(sf)
library(ozmaps)
source("constants.R")
library(ggrepel)
library(tidyverse)

load("cache/sites_data.rdata")
load("cache/metadata.rdata")
sea_dist <- read_rds("cache/sea_dist.rds")

maggie_kmls <- c("Huntingfield_Bay.kml","Wilson_Bay.kml","Maud_Bay.kml","Horseshoe_Bay.kml","Geoffrey_Bay.kml","Picnic_Bay.kml","Middle_Reef.kml")
extra_kmls <- c("Arthur_Bay.kml","Bay_Rock.kml","Florence_Bay.kml","Nelly_Bay.kml")

kml_files <- paste("dispersal",maggie_kmls,sep="/")
extra_kml_files <- paste("dispersal",extra_kmls,sep="/")

disp_data <- map(kml_files,st_read) %>% bind_rows()
disp_sf_sites <- disp_data %>% 
  mutate(name = str_replace(Name,"_"," ")) %>% 
  left_join(sites_data,by="name")

extra_sites <- readxl::read_excel("dispersal/extra_sites.xlsx") %>% 
  mutate(pop=NA,n_samples=NA,pop_order=NA,pop_group="None")



extra_disp_data <- map(extra_kml_files,st_read) %>% bind_rows()
extra_disp_sf_sites <- extra_disp_data %>% 
  mutate(name = str_replace(Name,"_"," ")) %>% 
  left_join(extra_sites, by="name")

gbr_coast <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "GBR_coastlineONLY")
gbr_features <- st_read("maps/GBR_NESP-TWQ-5/BaseLayers/",layer = "TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")

sf_use_s2(FALSE) # the input geometries are bad, and wrong spherically
x1 <- st_make_valid(gbr_features)

xmn=146.68
ymn=-19.25
xmx=146.9
ymx=-18.8

gbr_ft_cropped <- st_crop(x1, xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx) #%>% filter(LEVEL_1=="Island")
gbr_coast_cropped <- st_crop(gbr_coast,xmin = xmn, ymin = ymn, xmax = xmx, ymax = ymx)

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sites_sf <- st_as_sf(sites_data %>% filter(pop_group=="MI"),coords = c("lon", "lat"),crs = projcrs) %>% 
  mutate(name_label = paste(name," (",pop,")",sep=""))

extra_sites_sf <- st_as_sf(extra_sites,coords=c("lon","lat"),crs=projcrs) %>% 
  mutate(name_label = name)

sf_oz <- ozmap("states")

pmap <- ggplot() + 
  geom_sf(data = gbr_ft_cropped %>% filter(LEVEL_1=="Island")) +
  geom_sf(data = gbr_coast_cropped) +
  geom_sf(data = disp_sf_sites,aes(color=pop)) +
  geom_sf(data = extra_disp_sf_sites,color="grey") +
  geom_sf(data = sites_sf,aes(color=pop),size=3) + 
  geom_sf(data = extra_sites_sf,color="grey",size=3) +
  scale_color_manual(values = location_colors) +
  ggrepel::geom_text_repel(data = rbind(extra_sites_sf,sites_sf), mapping = aes(label=name_label, geometry=geometry), stat = "sf_coordinates",size=3) + 
#  ggrepel::geom_text_repel(data = extra_sites_sf, mapping = aes(label=name_label, geometry=geometry), stat = "sf_coordinates",size=3) +   
  theme_bw() +
  theme(legend.position="none") +
  xlab("") + ylab("") 

ggsave("figures/figure2A.svg",pmap,width = 120,units = "mm")

## IBD

ibd_plot <- sea_dist %>% 
  filter(comp_type!="between") %>% 
  filter(fst_type=="within_maggie") %>% 
  ggplot(aes(x=Dist,y=Fst)) + 
  geom_point(aes(color=finescale_pair_type)) + 
  geom_smooth(method="lm") +
  scale_color_manual(values = c("no_no"="grey50","so_no"="red","so_so"="grey20"), labels = c("no_no"="North vs North","so_so"="South vs South","so_no"="South vs North")) +
  theme_bw() +
  labs(
    color=NULL,
    x = "Sea Distance / km",
    y = expression("F"["st"])) +
  theme(legend.position = "bottom")

ggsave("figures/figure2B.svg",ibd_plot,width = 80,height=80,units = "mm")

# ## PCA
# 
load("cache/pcoa_ak_filtered_nr_mag.rdata")

pcoa_var_mag <- (pcoa_ak_filtered_nr_mag$eig/sum(pcoa_ak_filtered_nr_mag$eig))[1:2]

pcoa_ak_filtered_nr_meta <-  pcoa_ak_filtered_nr_mag$scores %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  left_join(metadata,by="ID") 

north <- c("HFB","WB","MB","HB")

south <- c("GB","PB","MR")

pca_plot <- pcoa_ak_filtered_nr_meta %>%
  mutate(side = case_when(
    pop %in% north ~ "North",
    pop %in% south ~ "South"
  )) %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=side)) +
  xlab(paste("PC1 (",round(pcoa_var_mag[1]*100,1)," %)",sep="")) +
  ylab(paste("PC2 (",round(pcoa_var_mag[2]*100,1)," %)",sep="")) + 
  theme_bw() +
  labs(color=NULL) +
  theme(legend.position = "bottom")

ggsave("figures/figure2C.svg",pca_plot,width = 80,height=80,units = "mm")

library(cowplot)

right_col <- plot_grid(
  ibd_plot+ theme(legend.position="none"),
  pca_plot+ theme(legend.position="none"),ncol=1)

bottom_row <- plot_grid(
  ibd_plot+ theme(legend.position="none"),
  pca_plot+ theme(legend.position="none"),ncol=2)

plot_grid(
  pmap + theme(legend.position="none"),
  bottom_row,
  ncol=1,
  rel_heights = c(4, 1))

plot_grid(
  pmap + theme(legend.position="none"),
  right_col,
  ncol=2,
  rel_widths = c(4, 1.5
                 ))

ggsave(filename = "figures/Figure2.png",width = 174,height = 110, units = "mm")
#ggsave(filename = "figures/Figure2.eps",width = 10,height = 8)
ggsave(filename = "figures/Figure2.svg",width = 174,height = 110, units = "mm")

