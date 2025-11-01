
# Named list to convert abbreviated names to long names

location_codes <- c("BRR","EP","GB","HaR","HB","HFB","JB","KR","LPB","MB","MR","PB","SO","WB","WP")
location_names <- c("Bramble Reef","East Pelorus","Geoffrey Bay","Havannah Reef","Horseshoe Bay","Huntingfield Bay","John Brewer Reef","Keeper Reef","Little Pioneer Bay","Maud Bay","Middle Reef","Picnic Bay","South Orpheus","Wilson Bay","West Pelorus")

names(location_names) <- location_codes

location_order <- 1:15
names(location_order) <- c("BRR","EP","WP","LPB","SO","JB","KR","HaR","HFB","WB","MB","HB","GB","PB","MR")

# Assign a certain colour to the clusters, so that cluster1 is consistently blue and cluster2 orange
cluster_colors <- c("blue","orange")

maggie_sites <- c("HFB","WB","MB","HB","GB","PB","MR")
maggie_no_sites <- c("HFB","WB","MB","HB")
maggie_so_sites <- c("GB","PB","MR")

## Assign colours to populations

#Assign colours based on these site locations: Magnetic Island vs. Palm Island vs. Mid-shelf

library(colorspace)

gr_or <- diverge_hcl(palette = "Green-Orange",n=8)

pop_colors = c("ADJ" = gr_or[1],"MI"=gr_or[8],"NO"=gr_or[7],"SO"=gr_or[6])

pop_names <- c("ADJ" = "Adjacent Reefs",
               "NO" = "North Magnetic Island",
               "SO" = "South Magnetic Island",
               "MI" = "Magnetic Island")
