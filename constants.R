
# Named list to convert abbreviated names to long names

location_codes <- c("BRR","EP","GB","HaR","HB","HFB","JB","KR","LPB","MB","MR","PB","SO","WB","WP")
location_names <- c("Bramble Reef","East Pelorus","Geoffrey Bay","Havannah Reef","Horseshoe Bay","Huntingfield Bay","John Brewer Reef","Keeper Reef","Little Pioneer Bay","Maud Bay","Middle Reef","Picnic Bay","South Orpheus","Wilson Bay","West Pelorus")

names(location_names) <- location_codes

location_order <- 1:15
names(location_order) <- c("BRR","EP","WP","LPB","SO","JB","KR","HaR","HFB","WB","MB","HB","GB","PB","MR")

# Assign a certain colour to the clusters, so that cluster1 is consistently blue and cluster2 orange
cluster_colors <- c("blue","orange")

maggie_sites <- c("HFB","WB","MB","HB","GB","PB","MR")

## Assign colours to populations

#Assign colours based on these site locations: Magnetic Island vs. Palm Island vs. Mid-shelf

# Get order of population names
#levels(ak.gen.mi@pop)
#  "BRR" "EP"  "MI"  "HaR" "JB"  "KR"  "LPB" "SO"  "WP"

# Assign colours to populations in that order:
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi <- c("#FDD835", "#74c476", "orange", "#74c476", "#FDD835", "#FDD835", "#74c476", "#74c476", "#74c476")

# Get order of population names
#levels(ak.gen.mi.pi@pop)
#  "BRR" "PI"  "MI"  "JB"  "KR"

# Assign colours to populations in that order: 
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi.pi <- c("#FDD835", "#74c476", "orange", "#FDD835", "#FDD835")

# Get order of population names
#levels(ak.gen.mi.pi.ms@pop)
#  "MS" "PI" "MI"

# Assign colours to populations in that order: 
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi.pi.ms <- c("#FDD835", "#74c476", "orange")

