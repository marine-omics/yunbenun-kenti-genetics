
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
